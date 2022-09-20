#include <iostream>
#include <iomanip>
#include <unordered_set>
#include <chrono>
#include <random>
#include <queue>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <xpas/phylo_kmer_db.h>
#include <xpas/version.h>
#include <xpas/phylo_tree.h>
#include <xpas/newick.h>
#include "db_builder.h"
#include "extended_tree.h"
#include "proba_matrix.h"
#include "ar.h"
#include "filter.h"
#include "branch_group.h"

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
        friend phylo_kmer_db build(const string& working_directory,
                                   const phylo_tree& original_tree, const phylo_tree& extended_tree,
                                   const proba_matrix& matrix,
                                   const ghost_mapping& mapping, const ar::mapping& ar_mapping,
                                   bool merge_branches, size_t kmer_size,
                                   phylo_kmer::score_type omega,
                                   filter_type filter, double mu, size_t num_threads);
    public:
        /// Member types


        /// \brief A group of node ids that must be processed together. We group together
        ///        the extended node ids that correspond to the same original node ids
        using id_group = std::vector<std::string>;

        /// \brief A group of probability submatrices that correspond to a group of nodes
        using proba_group = std::vector<std::reference_wrapper<const proba_matrix::mapped_type>>;


        /// Ctors, dtor and operator=
        db_builder(const string& working_directory,
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

        /// \brief The first stage of the construction algorithm. Creates a hashmap of phylo-kmers
        ///        for every node group.
        /// \return group ids (correspond to the post-order ids in the tree),
        ///         number of tuples explored, elapsed time
        std::tuple<std::vector<phylo_kmer::branch_type>, size_t, unsigned long> construct_group_hashmaps();

        /// \brief The second stage of the construction algorithm. Combines group hashmaps
        /// \return Elapsed time
        unsigned long merge_filtered(const std::vector<phylo_kmer::branch_type>& group_ids);

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
        std::pair<std::vector<group_hash_map>, size_t> explore_group(const proba_group& group) const;

        /// \brief Working and output directory
        string _working_directory;

        const phylo_tree& _original_tree;
        const phylo_tree& _extended_tree;

        const proba_matrix& _matrix;
        const ghost_mapping& _extended_mapping;
        const ar::mapping& _ar_mapping;

        bool _merge_branches;

        size_t _kmer_size;
        xpas::phylo_kmer::score_type _omega;

        double _reduction_ratio;

        xpas::filter_type _filter;
        double _mu;

        /// The number of batches in which the space of k-mers is split
        const size_t _num_batches = 1;

        size_t _num_threads;
        phylo_kmer_db _phylo_kmer_db;

    };

    db_builder::db_builder(const string& working_directory,
                           const phylo_tree& original_tree, const phylo_tree& extended_tree,
                           const proba_matrix& matrix,
                           const ghost_mapping& mapping, const ar::mapping& ar_mapping,
                           bool merge_branches, size_t kmer_size, phylo_kmer::score_type omega,
                           filter_type filter, double mu, size_t num_threads)
        : _working_directory{ working_directory }
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
        std::cout << "Construction parameters:" << std::endl <<
                  "\tSequence type: " << xpas::seq_type::name << std::endl <<
                  "\tk: " << _kmer_size << std::endl <<
                  "\tomega: " << _omega << std::endl <<
                  "\tKeep positions: " << (xpas::keep_positions ? "true" : "false") << std::endl << std::endl;

        /// The first stage of the algorithm: create a hashmap for every node group
        std::cout << "Building database: [stage 1 / 2]:" << std::endl;
        const auto& [group_ids, num_tuples, construction_time] = construct_group_hashmaps();
        std::cout << "Calculated " << num_tuples << " phylo-k-mers.\nCalculation time: " << construction_time
                  << "\n\n" << std::flush;

        /// The second stage of the algorithm: combine merge and filter
        std::cout << "Building database: [stage 2 / 2]:" << std::endl;
        const auto merge_time = merge_filtered(group_ids);
        fs::remove_all(get_groups_dir(_working_directory));

        /// Calculate the number of phylo-kmers stored in the database
        size_t total_entries = 0;
        for (const auto& kmer_entry : _phylo_kmer_db)
        {
            total_entries += kmer_entry.second.size();
        }

        std::cout << "Building database: Done.\n";
        std::cout << "Built " << total_entries << " phylo-k-mers for "
                  << _phylo_kmer_db.size() << " different k-mers.\nTotal time (ms): "
                  << construction_time + merge_time << "\n\n" << std::flush;
    }

    std::tuple<std::vector<phylo_kmer::branch_type>, size_t, unsigned long> db_builder::construct_group_hashmaps()
    {
        /// create a temporary directory for hashmaps
        const auto temp_dir = get_groups_dir(_working_directory);
        fs::create_directories(temp_dir);

        try
        {
            /// Run the construction algorithm
            const auto begin = std::chrono::steady_clock::now();
            const auto& [group_ids, num_tuples] = explore_kmers();
            const auto end = std::chrono::steady_clock::now();
            const auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
            return { group_ids, num_tuples, elapsed_time };
        }
        catch (const std::exception& error)
        {
            std::cerr << "Error: " << error.what() << std::endl;
            fs::remove_all(temp_dir);
            throw error;
        }
    }


    unsigned long db_builder::merge_filtered(const std::vector<phylo_kmer::branch_type>& group_ids)
    {
        const auto begin = std::chrono::steady_clock::now();

        /// Filter phylo k-mers
        const auto threshold = xpas::score_threshold(_omega, _kmer_size);
        auto filter = xpas::make_filter(_filter, _original_tree.get_node_count(),
                                        _working_directory, _num_batches, _mu, threshold);
        filter->filter(group_ids);

        size_t filtered_kmers = 0;
        size_t total_kmers = 0;
        size_t total_entries = 0;
        size_t filtered_entries = 0;
        for (size_t batch_idx = 0; batch_idx < _num_batches; ++batch_idx)
        {
            const auto temp_db = xpas::merge_batch(_working_directory, group_ids, batch_idx);
            for (const auto& [key, entries] : temp_db)
            {
                total_entries += entries.size();
                total_kmers++;
                if (filter->is_good(key))
                {
                    filtered_kmers++;
                    filtered_entries += entries.size();
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
                        //std::cout << key << " " << xpas::decode_kmer(key, _kmer_size) << ": " << std::endl;
                        for (const auto& [branch, score] : entries)
                        {
                            //std::cout << "\t\t" << branch << " -> " << score << " " << std::pow(10, score) << std::endl;
                            _phylo_kmer_db.unsafe_insert(key, {branch, score});
                        }
                    }
#endif
                }
            }
        }

        const auto end = std::chrono::steady_clock::now();
        const auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

        std::cout << "Kept " << filtered_kmers << " / " << total_kmers
                  << " k-mers (" << std::setprecision(3) << ((float) filtered_kmers) / total_kmers * 100 << "%) | "
                  << filtered_entries << " / " << total_entries
                  << " entries (" << std::setprecision(3) << ((float) filtered_entries) / total_entries * 100 << "%)."
                  << "\nFiltering time: " << time << "\n\n" << std::flush;
        return time;
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


            /// Get sub-matrices of probabilities for the group
            const auto matrices = get_submatrices(node_group);

            /// Explore k-mers of the group and store results in a hash map
            const auto& [hash_maps, branch_count] = explore_group(matrices);

            /// Save the group hashmap on disk
            size_t index = 0;
            for (const auto& hash_map : hash_maps)
            {
                save_group_map(hash_map, get_group_map_file(_working_directory, original_node_postorder_id, index));
                ++index;
            }

            count += branch_count;
        }
        return { node_postorder_ids, count };
    }

    #ifdef KEEP_POSITIONS
    std::pair<std::vector<group_hash_map>, size_t> db_builder::explore_group(const proba_group& group) const
    {
        auto hash_maps = std::vector<group_hash_map>(_num_batches);
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
                    put(hash_maps[kmer_batch(kmer.key, _num_batches)], positioned_kmer);
                    ++count;
                }
            }
        }

        return {std::move(hash_maps), count};
    }
    #else

    std::pair<std::vector<group_hash_map>, size_t> db_builder::explore_group(const proba_group& group) const
    {
        auto hash_maps = std::vector<group_hash_map>(_num_batches);

        size_t count = 0;
        const auto log_threshold = std::log10(xpas::score_threshold(_omega, _kmer_size));
        for (auto node_entry_ref : group)
        {
            const auto& node_entry = node_entry_ref.get();

            for (auto& window : chain_windows(node_entry, _kmer_size, log_threshold))
            {
                for (const auto& kmer : window)
                {
                    xpas::put(hash_maps[kmer_batch(kmer.key, _num_batches)], kmer);
                    ++count;
                }
            }
        }

        return {std::move(hash_maps), count};
    }

    #endif
};


namespace xpas
{

    phylo_kmer_db build(const string& working_directory,
                        const phylo_tree& original_tree, const phylo_tree& extended_tree,
                        const proba_matrix& matrix,
                        const ghost_mapping& mapping, const ar::mapping& ar_mapping,
                        bool merge_branches,
                        size_t kmer_size, xpas::phylo_kmer::score_type omega, filter_type filter, double mu, size_t num_threads)
    {
        db_builder builder(working_directory,
                           original_tree, extended_tree,
                           matrix,
                           mapping, ar_mapping,
                           merge_branches, kmer_size, omega, filter, mu, num_threads);
        builder.run();
        return std::move(builder._phylo_kmer_db);
    }
}
