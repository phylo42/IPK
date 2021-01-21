#include <iostream>
#include <unordered_set>
#include <chrono>
#include <random>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <xpas/phylo_kmer_db.h>
#include <xpas/version.h>
#include <xpas/phylo_tree.h>
#include <xpas/newick.h>
#include "db_builder.h"
#include "extended_tree.h"
#include "alignment.h"
#include "proba_matrix.h"
#include "ar.h"

using std::string;
using std::cout, std::endl;
using std::to_string;
using namespace xpas;
namespace fs = boost::filesystem;


namespace xpas
{
    /// \brief Constructs a database of phylo-kmers.
    class db_builder
    {
        friend phylo_kmer_db build(string working_directory,
                                   const alignment& original_alignment, const alignment& extended_alignment,
                                   const phylo_tree& original_tree, const phylo_tree& extended_tree,
                                   const proba_matrix& matrix,
                                   const ghost_mapping& mapping, const ar::mapping& ar_mapping,
                                   bool merge_branches, size_t kmer_size,
                                   phylo_kmer::score_type omega,
                                   filter_type filter, double mu, size_t num_threads);
    public:
        /// Member types
        /// \brief A hash map to store all the phylo-kmers, placed to one original node
#ifdef KEEP_POSITIONS

        struct score_pos_pair
        {
            phylo_kmer::score_type score;
            phylo_kmer::pos_type position;
        };

        using branch_hash_map = hash_map<phylo_kmer::key_type, score_pos_pair>;
#else
        using branch_hash_map = hash_map<phylo_kmer::key_type, phylo_kmer::score_type>;
#endif

        /// \brief A group of node ids that must be processed together. We group together
        ///        the extended node ids that correspond to the same original node ids
        using id_group = std::vector<std::string>;

        /// \brief A group of probability submatrices that correspond to a group of nodes
        using proba_group = std::vector<std::reference_wrapper<const proba_matrix::mapped_type>>;


        /// Ctors, dtor and operator=
        db_builder(string working_directory,
                   const alignment& original_alignment, const alignment& extended_alignment,
                   const phylo_tree& original_tree, const phylo_tree& extended_tree,
                   const proba_matrix& matrix,
                   const ghost_mapping& mapping, const ar::mapping& ar_mapping,
                   bool merge_branches, size_t kmer_size, phylo_kmer::score_type omega,
                   filter_type filter, double mu, size_t num_threads);
        db_builder(const db_builder&) = delete;
        db_builder(db_builder&&) = delete;
        db_builder& operator=(const db_builder&) = delete;
        db_builder& operator=(db_builder&&) = delete;
        ~db_builder() noexcept = default;


        /// \brief Runs the database construction
        void run();

    private:
        /// \brief Results of the first stage of the algorithm:
        using explore_groups_result = std::tuple<std::vector<phylo_kmer::branch_type>, size_t, unsigned long>;

        /// \brief The first stage of the construction algorithm. Creates a hashmap of phylo-kmers
        ///        for every node group.
        /// \return group ids (correspond to the post-order ids in the tree),
        ///         number of tuples explored, elapsed time
        std::tuple<std::vector<phylo_kmer::branch_type>, size_t, unsigned long> construct_group_hashmaps();

        /// \brief The second stage of the construction algorithm. Combines group hashmaps
        /// \return Elapsed time
        unsigned long merge_and_filter(const std::vector<phylo_kmer::branch_type>& group_ids);

        /// Merges hashmaps of the same index into a database
        [[nodiscard]]
        phylo_kmer_db merge_hashmaps(const std::vector<phylo_kmer::branch_type>& group_ids, size_t index) const;

        [[nodiscard]]
        hash_map<phylo_kmer::key_type, bool> filter_keys(const phylo_kmer_db& db) const;

        /// \brief Returns a filename for a hashmap of a given group
        [[nodiscard]]
        std::string group_hashmap_file(const branch_type& group, size_t index) const;

        /// \brief Groups ghost nodes by corresponding original node id
        [[nodiscard]]
        std::vector<id_group> group_ghost_ids(const std::vector<std::string>& ghost_ids) const;

        /// \brief Groups references to submatrices of probabilites, corresponding to a group of nodes
        [[nodiscard]]
        proba_group get_submatrices(const id_group& group) const;

        /// \brief Runs a phylo-kmer exploration for every ghost node of the extended_tree
        /// \return 1) A vector of group ids, which correspond to post-order node ids in the tree
        ///         2) The number of explored phylo-kmers. This number can be more than a size of a resulting database
        [[nodiscard]]
        std::tuple<std::vector<phylo_kmer::branch_type>, size_t> explore_kmers() const;

        /// \brief Explores phylo-kmers of a collection of ghost nodes. Here we assume that the nodes
        ///        in the group correspond to one original node
        /// \return A hash map with phylo-kmers stored and a number of explored phylo-kmers
        [[nodiscard]]
        std::pair<std::vector<branch_hash_map>, size_t> explore_group(const proba_group& group) const;

        /// \brief Saves a hash map to file
        void save_hash_map(const branch_hash_map& map, const std::string& filename) const;

        /// \brief Loads a hash map from file
        [[nodiscard]]
        branch_hash_map load_hash_map(const std::string& filename) const;

        /// \brief Working and output directory
        string _working_directory;
        /// \brief A subdirectory of _working_directory to store hashmaps
        string _hashmaps_directory;

        const alignment& _original_alignment;
        const alignment& _extended_alignment;

        const phylo_tree& _original_tree;
        const phylo_tree& _extended_tree;

        const proba_matrix& _matrix;
        const ghost_mapping& _extended_mapping;
        const ar::mapping& _ar_mapping;

        bool _merge_branches;

        size_t _kmer_size;
        xpas::phylo_kmer::score_type _omega;

        double _reduction_ratio;

        filter_type _filter;
        double _mu;

        /// The number of ranges in which the space of k-mers is split
        const size_t _num_ranges = 4;

        size_t _num_threads;
        phylo_kmer_db _phylo_kmer_db;

    };

    db_builder::db_builder(string working_directory,
                           const alignment& original_alignment, const alignment& extended_alignment,
                           const phylo_tree& original_tree, const phylo_tree& extended_tree,
                           const proba_matrix& matrix,
                           const ghost_mapping& mapping, const ar::mapping& ar_mapping,
                           bool merge_branches, size_t kmer_size, phylo_kmer::score_type omega,
                           filter_type filter, double mu, size_t num_threads)
        : _working_directory{ std::move(working_directory) }
        , _hashmaps_directory{ (fs::path{working_directory} / fs::path{"hashmaps"}).string() }
        , _original_alignment{ original_alignment }
        , _extended_alignment{ extended_alignment }
        , _original_tree{ original_tree }
        , _extended_tree{ extended_tree }
        , _matrix{ matrix }
        , _extended_mapping{ mapping }
        , _ar_mapping{ ar_mapping }
        , _merge_branches{ merge_branches }
        , _kmer_size{ kmer_size }
        , _omega{ omega }
        , _filter{ filter }
        , _mu{ mu }
        , _num_threads{ num_threads }
        , _phylo_kmer_db{ kmer_size, omega, xpas::seq_type::name, xpas::io::to_newick(_original_tree)}
    {}

    void db_builder::run()
    {
        /// Run ancestral reconstruction

        /// The first stage of the algorithm -- create a hashmap for every group node
        const auto& [group_ids, num_tuples, construction_time] = construct_group_hashmaps();

        /// The second stage of the algorithm -- combine hashmaps
        const auto merge_time = merge_and_filter(group_ids);

        /// Calculate the number of phylo-kmers stored in the database
        size_t total_entries = 0;
        for (const auto& kmer_entry : _phylo_kmer_db)
        {
            total_entries += kmer_entry.second.size();
        }

        std::cout << "Built " << total_entries << " phylo-kmers out of " << num_tuples << " for "
                  << _phylo_kmer_db.size() << " k-mer values.\nTime (ms): "
                  << construction_time + merge_time << "\n\n" << std::flush;
    }

    std::tuple<std::vector<phylo_kmer::branch_type>, size_t, unsigned long> db_builder::construct_group_hashmaps()
    {
        /// create a temporary directory for hashmaps
        fs::create_directories(_hashmaps_directory);

        std::cout << "Construction parameters:" << std::endl <<
                     "\tSequence type: " << xpas::seq_type::name << std::endl <<
                     "\tk: " << _kmer_size << std::endl <<
                     "\tomega: " << _omega << std::endl <<
                     "\tKeep positions: " << (xpas::keep_positions ? "true" : "false") << std::endl << std::endl;

        std::cout << "Building database..." << std::endl;

        /// Run the construction algorithm
        const auto begin = std::chrono::steady_clock::now();
        const auto& [group_ids, num_tuples] = explore_kmers();
        const auto end = std::chrono::steady_clock::now();
        const auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        return { group_ids, num_tuples, elapsed_time };
    }

    double shannon(double x)
    {
        return - x * std::log2(x);
    }

    xpas::phylo_kmer::score_type logscore_to_score(xpas::phylo_kmer::score_type log_score)
    {
        return std::min(std::pow(10, log_score), 1.0);
    }

    unsigned long db_builder::merge_and_filter(const std::vector<phylo_kmer::branch_type>& group_ids)
    {
        const auto begin = std::chrono::steady_clock::now();

        if (_filter == filter_type::entropy)
        {
            std::cout << "Filtering (minimal conditional entropy)..." << std::endl;
        }
        else if (_filter == filter_type::random)
        {
            std::cout << "Filtering (random)..." << std::endl;
        }
        else
        {
            std::cout << "No filtering." << std::endl;
        }

        for (size_t index = 0; index < _num_ranges; ++index)
        {
            auto temp_db = merge_hashmaps(group_ids, index);
            hash_map<phylo_kmer::key_type, bool> keep_keys = filter_keys(temp_db);

            for (const auto& [key, entries] : temp_db)
            {
                if (keep_keys[key])
                {
#ifdef KEEP_POSITIONS
                    if (_merge_branches)
                    {
                        for (const auto& [branch, score, position] : entries)
                        {
                            if (auto entries = _phylo_kmer_db.search(key); entries)
                            {
                                /// If there are entries, there must be only one because
                                /// we always take maximum score among different branches.
                                /// So this loop will have only one iteration
                                for (const auto& [old_node, old_score, old_position] : *entries)
                                {
                                    if (old_score < score)
                                    {
                                        _phylo_kmer_db.replace(key, { branch, score, position });
                                    }
                                }
                            }
                            else
                            {
                                _phylo_kmer_db.unsafe_insert(key, { branch, score, position });
                            }
                        }
                    }
                    else
                    {
                        for (const auto& [branch, score, position] : entries)
                        {
                            _phylo_kmer_db.unsafe_insert(key, { branch, score, position });
                        }
                    }
#else
                    if (_merge_branches)
                    {
                        throw std::runtime_error("--merge-branches is only supported for xpas compiled with the KEEP_POSITIONS flag.");
                    }
                    else
                    {
                        for (const auto& [branch, score] : entries)
                        {
                            _phylo_kmer_db.unsafe_insert(key, {branch, score});
                        }
                    }
#endif
                }
            }

            size_t counter = 0;
            for (const auto& [_, keep] : keep_keys)
            {
                if (keep)
                {
                    counter++;
                }
            }
            std::cout << counter << " out of " << keep_keys.size() <<
                      " (" << ((float) counter) / keep_keys.size() << ")" << std::endl;
        }


        const auto end = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    }

    bool is_ghost(const phylo_node& node)
    {
        const string& label = node.get_label();
        return boost::ends_with(label, "_X0") || boost::ends_with(label, "_X1");
    }

    /// \brief Returns a list of ghost node ids
    std::vector<std::string> get_ghost_ids(const phylo_tree& tree)
    {
        std::vector<std::string> branch_ids;

        for (const auto& branch_node: tree)
        {
            if (is_ghost(branch_node))
            {
                branch_ids.push_back(branch_node.get_label());
            }
        }
        return branch_ids;
    }

    phylo_kmer_db db_builder::merge_hashmaps(const std::vector<phylo_kmer::branch_type>& group_ids, size_t index) const
    {
        std::cout << "Merging hash maps [index = " << index << "]..." << std::endl;
        phylo_kmer_db temp_db(_phylo_kmer_db.kmer_size(), _phylo_kmer_db.omega(), xpas::seq_type::name, "");

        /// Load hash maps and merge them
        for (const auto group_id : group_ids)
        {
            const auto hash_map = load_hash_map(group_hashmap_file(group_id, index));
#ifdef KEEP_POSITIONS
            for (const auto& [key, score_pos_pair] : hash_map)
            {
                const auto& [score, position] = score_pos_pair;
                temp_db.unsafe_insert(key, {group_id, score, position});
            }
#else
            for (const auto& [key, score] : hash_map)
            {
                temp_db.unsafe_insert(key, {group_id, score});
            }
#endif
        }

        return temp_db;
    }

    hash_map<phylo_kmer::key_type, bool> db_builder::filter_keys(const phylo_kmer_db& db) const
    {
        std::cout << "Filtering phylo k-mers..." << std::endl;

        /// create a hash map which contains boolean values
        /// which tell if we keep corresponding k-mers.
        hash_map<phylo_kmer::key_type, bool> keep_keys;
        for (const auto& [key, _] : db)
        {
            keep_keys[key] = true;
        }

        if (_filter == filter_type::entropy)
        {
            const auto threshold = xpas::score_threshold(_phylo_kmer_db.omega(), _phylo_kmer_db.kmer_size());
            const auto log_threshold = std::log10(threshold);

            hash_map<phylo_kmer::key_type, double> filter_stats;
            for (const auto& [key, entries] : db)
            {
                /// calculate the score sum to normalize scores
                double score_sum = 0;
                double log_score_sum = 0;
#ifdef KEEP_POSITIONS
                for (const auto& [branch, log_score, position] : entries)
                {
                    (void)branch;
                    (void)position;
                    score_sum += logscore_to_score(log_score);
                    log_score_sum += log_score;
                }
#else
                for (const auto& [_, log_score] : entries)
                {
                    score_sum += logscore_to_score(log_score);
                    log_score_sum += log_score;
                }
#endif

                /// do not forget the branches that are not stored in the database,
                /// they suppose to have the threshold score
                score_sum += static_cast<double>(_original_tree.get_node_count() - entries.size()) * threshold;
                log_score_sum +=
                    static_cast<double>(_original_tree.get_node_count() - entries.size()) * log_threshold;

                /// Entropy
                const auto weighted_threshold = threshold / score_sum;
                const auto target_threshold = shannon(weighted_threshold);

#ifdef KEEP_POSITIONS
                for (const auto& [branch, log_score, position] : entries)
#else
                for (const auto& [branch, log_score] : entries)
#endif
                {
                    const auto weighted_score = logscore_to_score(log_score) / score_sum;
                    const auto target_value = shannon(weighted_score);

                    if (filter_stats.find(key) == filter_stats.end())
                    {
                        filter_stats[key] = static_cast<double>(_original_tree.get_node_count()) * target_threshold;
                    }
                    filter_stats[key] = filter_stats[key] - target_threshold + target_value;
                }
            }

            /*for (const auto& [key, stats] : filter_stats)
            {
                std::cout << stats << " ";
            }
            std::cout << std::endl << std::endl;
*/

            /// copy entropy values into a vector
            std::vector<phylo_kmer::score_type> filter_values;
            filter_values.reserve(filter_stats.size());
            for (const auto& entry : filter_stats)
            {
                filter_values.push_back(entry.second);
            }

            std::cout << "Partial sort..." << std::endl;
            const auto qth_element = size_t(filter_values.size() * _mu);
            std::nth_element(filter_values.begin(), filter_values.begin() + qth_element, filter_values.end());
            const auto quantile = filter_values[qth_element];

            std::cout << "Filter threshold value: " << quantile << std::endl;

            for (const auto& [key, entries] : db)
            {
                keep_keys[key] = (filter_stats[key] <= quantile);
            }
        }
        else if (_filter == filter_type::random)
        {
            std::default_random_engine generator;
            std::uniform_real_distribution<double> distribution(0, 1);

            /// update keep_keys randomly
            for (const auto& [key, _] : db)
            {
                keep_keys[key] = distribution(generator) <= _mu;
            }
        }
        else
        {
            ///  No filtering
        }

        return keep_keys;
    }

    std::string db_builder::group_hashmap_file(const branch_type& group, size_t index) const
    {
        return (fs::path{_hashmaps_directory} / fs::path{std::to_string(group) + "_" + std::to_string(index)}).string();
    }

    std::vector<db_builder::id_group> db_builder::group_ghost_ids(const std::vector<std::string>& ghost_ids) const
    {
        std::vector<id_group> groups;
        groups.reserve(ghost_ids.size() / 2);

        std::unordered_map<branch_type, size_t> index_mapping;
        for (const auto& ghost_id : ghost_ids)
        {
            const auto original_postorder_id = _extended_mapping.at(ghost_id);

            /// Ignore the root
            const auto original_node = _original_tree.get_by_postorder_id(original_postorder_id);
            if (original_node && (*original_node)->is_root())
            {
                continue;
            }

            if (const auto it = index_mapping.find(original_postorder_id); it != index_mapping.end())
            {
                groups[it->second].push_back(ghost_id);
            }
            else
            {
                groups.push_back({ghost_id});
                index_mapping[original_postorder_id] = groups.size() - 1;
            }

        }
        return groups;
    }

    db_builder::proba_group db_builder::get_submatrices(const id_group& group) const
    {
        proba_group submatrices;

        for (const auto& ext_node_label : group)
        {
            const auto& ar_node_label = _ar_mapping.at(ext_node_label);
            if (const auto& it = _matrix.find(ar_node_label); it != _matrix.end())
            {
                submatrices.push_back(std::cref(it->second));
            }
            else
            {
                throw std::runtime_error("Internal error: could not find " + ar_node_label + " node. "
                                         "Make sure it is in the ARTree_id_mapping file.");
            }

        }

        return submatrices;
    }

    std::tuple<std::vector<phylo_kmer::branch_type>, size_t> db_builder::explore_kmers() const
    {
        size_t count = 0;

        /// Here we assume that every original node corresponds to two ghost nodes.
        const size_t ghosts_per_node = 2;

        /// Filter and group ghost nodes
        const auto node_groups = group_ghost_ids(get_ghost_ids(_extended_tree));

        /// Process branches in parallel. Results of the branch-and-bound algorithm are stored
        /// in a hash map for every group separately on disk.
        std::vector<phylo_kmer::branch_type> node_postorder_ids(node_groups.size());

        /*
        #pragma omp parallel for schedule(auto) reduction(+: count) num_threads(_num_threads) \
            shared(original_tree, node_postorder_ids, _matrix)
        */
        for (size_t i = 0; i < node_groups.size(); ++i)
        {
            const auto& node_group = node_groups[i];
            assert(node_group.size() == ghosts_per_node);
            (void) ghosts_per_node;

            /// Having a label of a node in the extended tree, we need to find the corresponding node
            /// in the original tree. We take the first ghost node, because all of them correspond to
            /// the same original node
            const auto original_node_postorder_id = _extended_mapping.at(node_group[0]);
            node_postorder_ids[i] = original_node_postorder_id;

            /// Get sub-matrices of probabilities for a group
            const auto matrices = get_submatrices(node_group);

            /// Explore k-mers of a group and store results in a hash map
            const auto& [hash_maps, branch_count] = explore_group(matrices);

            /// Save hash maps on disk
            size_t index = 0;
            for (const auto& hash_map : hash_maps)
            {
                save_hash_map(hash_map, group_hashmap_file(original_node_postorder_id, index));
                ++index;
            }

            count += branch_count;
        }
        return { node_postorder_ids, count };
    }

    /// \brief Puts a key-value pair in a hash map. Used to process branches in parallel
    #ifdef KEEP_POSITIONS
    void put(db_builder::branch_hash_map& map, const phylo_kmer& kmer)
    {
        if (auto it = map.find(kmer.key); it != map.end())
        {
            if (it->second.score < kmer.score)
            {
                map[kmer.key] = { kmer.score, kmer.position };
            }
        }
        else
        {
            map[kmer.key] = { kmer.score, kmer.position };
        }
    }

    std::pair<std::vector<db_builder::branch_hash_map>, size_t> db_builder::explore_group(const proba_group& group) const
    {
        auto hash_maps = std::vector<branch_hash_map>(_num_ranges);
        size_t count = 0;

        const auto log_threshold = std::log10(xpas::score_threshold(_omega, _kmer_size));
        for (auto node_entry_ref : group)
        {
            const auto& node_entry = node_entry_ref.get();
            for (auto& window : chain_windows(node_entry, _kmer_size, log_threshold))
            {
                const auto position = window.get_start_pos();
                for (const auto& kmer : window)
                {
                    phylo_kmer positioned_kmer = { kmer.key, kmer.score, position };
                    put(hash_maps[kmer.key % _num_ranges], positioned_kmer);
                    ++count;
                }
            }
        }

        return {std::move(hash_maps), count};
    }
    #else

    void put(db_builder::branch_hash_map& map, const phylo_kmer& kmer)
    {
        if (auto it = map.find(kmer.key); it != map.end())
        {
            if (it->second < kmer.score)
            {
                map[kmer.key] = kmer.score;
            }
        }
        else
        {
            map[kmer.key] = kmer.score;
        }
    }

    std::pair<std::vector<db_builder::branch_hash_map>, size_t> db_builder::explore_group(const proba_group& group) const
    {
        auto hash_maps = std::vector<branch_hash_map>(_num_ranges);

        size_t count = 0;
        const auto log_threshold = std::log10(xpas::score_threshold(_omega, _kmer_size));
        for (auto node_entry_ref : group)
        {
            const auto& node_entry = node_entry_ref.get();

            //std::cout << node_entry.get_label() << std::endl;
            for (auto& window : chain_windows(node_entry, _kmer_size, log_threshold))
            {
                //std::cout << "\tWINDOW " << window.get_start_pos() << std::endl;
                for (const auto& kmer : window)
                {
                    //std::cout << "\t\t" << kmer.key << " " << xpas::decode_kmer(kmer.key, _kmer_size) << " -> "
                    //          << kmer.score << " " << std::pow(10, kmer.score) << std::endl;
                    put(hash_maps[kmer.key % _num_ranges], kmer);
                    ++count;
                }
            }
        }

        return {std::move(hash_maps), count};
    }

    #endif
}

namespace boost::serialization
{

#ifdef KEEP_POSITIONS
    /// Serialize a hash map
    template<class Archive>
    inline void save(Archive& ar, const ::db_builder::branch_hash_map& map, const unsigned int /*version*/)
    {
        size_t map_size = map.size();
        ar & map_size;

        for (const auto& [key, score_pos_pair] : map)
        {
            const auto& [score, position] = score_pos_pair;
            ar & key & score & position;
        }
    }

    /// Deserialize a hash map
    template<class Archive>
    inline void load(Archive& ar, ::db_builder::branch_hash_map& map, unsigned int version)
    {
        (void)version;

        size_t map_size = 0;
        ar & map_size;

        for (size_t i = 0; i < map_size; ++i)
        {
            xpas::phylo_kmer::key_type key = xpas::phylo_kmer::na_key;
            xpas::phylo_kmer::score_type score = xpas::phylo_kmer::na_score;
            xpas::phylo_kmer::pos_type position = xpas::phylo_kmer::na_pos;
            ar & key & score & position;
            map[key] = { score, position };
        }
    }
#else
    /// Serialize a hash map
    template<class Archive>
    inline void save(Archive& ar, const ::db_builder::branch_hash_map& map, const unsigned int /*version*/)
    {
        size_t map_size = map.size();
        ar & map_size;

        for (const auto&[key, score] : map)
        {
            ar & key & score;
        }
    }

    /// Deserialize a hash map
    template<class Archive>
    inline void load(Archive& ar, ::db_builder::branch_hash_map& map, unsigned int version)
    {
        (void)version;

        size_t map_size = 0;
        ar & map_size;

        for (size_t i = 0; i < map_size; ++i)
        {
            xpas::phylo_kmer::key_type key = xpas::phylo_kmer::na_key;
            xpas::phylo_kmer::score_type score = xpas::phylo_kmer::na_score;
            ar & key & score;
            map[key] = score;
        }
    }
#endif

    // split non-intrusive serialization function member into separate
    // non intrusive save/load member functions
    template<class Archive>
    inline void serialize(Archive& ar, ::db_builder::branch_hash_map& map, unsigned int file_version)
    {
        boost::serialization::split_free(ar, map, file_version);
    }
}

namespace xpas
{
    void db_builder::save_hash_map(const branch_hash_map& map, const std::string& filename) const
    {
        std::ofstream ofs(filename);
        boost::archive::binary_oarchive oa(ofs);
        oa & map;
    }

    /// \brief Loads a hash map from file
    db_builder::branch_hash_map db_builder::load_hash_map(const std::string& filename) const
    {
        std::ifstream ifs(filename);
        boost::archive::binary_iarchive ia(ifs);

        ::db_builder::branch_hash_map map;
        ia & map;
        return map;
    }

    phylo_kmer_db build(string working_directory,
                        const alignment& original_alignment, const alignment& extended_alignment,
                        const phylo_tree& original_tree, const phylo_tree& extended_tree,
                        const proba_matrix& matrix,
                        const ghost_mapping& mapping, const ar::mapping& ar_mapping,
                        bool merge_branches,
                        size_t kmer_size, xpas::phylo_kmer::score_type omega, filter_type filter, double mu, size_t num_threads)
    {
        db_builder builder(std::move(working_directory),
                           original_alignment, extended_alignment,
                           original_tree, extended_tree,
                           matrix,
                           mapping, ar_mapping,
                           merge_branches, kmer_size, omega, filter, mu, num_threads);
        builder.run();
        return std::move(builder._phylo_kmer_db);
    }
}
