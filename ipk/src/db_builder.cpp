#include <iostream>
#include <iomanip>
#include <unordered_set>
#include <chrono>
#include <random>
#include <queue>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>
#include <i2l/phylo_kmer_db.h>
#include <i2l/serialization.h>
#include <i2l/version.h>
#include <i2l/phylo_tree.h>
#include <i2l/newick.h>
#include "db_builder.h"
#include "extended_tree.h"
#include "proba_matrix.h"
#include "ar.h"
#include "filter.h"
#include "branch_group.h"
#include "pk_compute.h"


using std::string;
using std::cout, std::endl;
using std::to_string;
using namespace i2l;
using namespace ipk;
namespace fs = boost::filesystem;


namespace ipk
{
    /// \brief Constructs a database of phylo-kmers.
    class db_builder
    {
        friend phylo_kmer_db build(const string& working_directory,
                                   const phylo_tree& original_tree, const phylo_tree& extended_tree,
                                   const proba_matrix& matrix,
                                   const ghost_mapping& mapping, const ar::mapping& ar_mapping, bool merge_branches,
                                   ipk::algorithm algorithm, ipk::ghost_strategy strategy,
                                   size_t kmer_size, phylo_kmer::score_type omega,
                                   filter_type filter, double mu,
                                   size_t num_threads);
    public:
        /// Member types


        /// \brief A group of node ids that must be processed together. We group together
        ///        the extended node ids that correspond to the same original node ids
        using id_group = std::vector<std::string>;

        /// \brief A group of probability submatrices that correspond to a group of nodes
        using proba_group = std::vector<std::reference_wrapper<const proba_matrix::mapped_type>>;


        /// Ctors, dtor and operator=
        db_builder(std::string working_directory,
                   const phylo_tree& original_tree, const phylo_tree& extended_tree,
                   const proba_matrix& matrix,
                   const ghost_mapping& mapping, const ar::mapping& ar_mapping, bool merge_branches,
                   ipk::algorithm algorithm, ipk::ghost_strategy strategy,
                   size_t kmer_size, phylo_kmer::score_type omega,
                   filter_type filter, double mu,
                   size_t num_threads);
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
        /// and runs batched filtering
        /// \return Elapsed time
        unsigned long merge_filtered(const std::vector<phylo_kmer::branch_type>& group_ids);

        /// \brief The second stage of the construction algorithm. Merges group hashmaps
        unsigned long merge(const std::vector<phylo_kmer::branch_type>& group_ids);

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
        std::pair<std::vector<group_hash_map>, size_t> explore_group(const proba_group& group, size_t postorder_id) const;

        /// \brief Working and output directory
        string _working_directory;

        const phylo_tree& _original_tree;
        const phylo_tree& _extended_tree;

        //const proba_matrix& _matrix;
        const proba_matrix& _matrix;
        const ghost_mapping& _extended_mapping;
        const ar::mapping& _ar_mapping;

        bool _merge_branches;

        ipk::algorithm _algorithm;
        ipk::ghost_strategy _ghost_strategy;

        size_t _kmer_size;
        i2l::phylo_kmer::score_type _omega;

        ipk::filter_type _filter;
        double _mu;

        /// The number of batches in which the space of k-mers is split
        const size_t _num_batches = 4;

        size_t _num_threads;
        phylo_kmer_db _phylo_kmer_db;
    };

    db_builder::db_builder(std::string working_directory,
                           const phylo_tree& original_tree, const phylo_tree& extended_tree,
                           const proba_matrix& matrix,
                           const ghost_mapping& mapping, const ar::mapping& ar_mapping, bool merge_branches,
                           ipk::algorithm algorithm, ipk::ghost_strategy strategy,
                           size_t kmer_size, phylo_kmer::score_type omega,
                           filter_type filter, double mu,
                           size_t num_threads)
        : _working_directory{ std::move(working_directory) }
        , _original_tree{ original_tree }
        , _extended_tree{ extended_tree }
        , _matrix{ matrix }
        , _extended_mapping{ mapping }
        , _ar_mapping{ ar_mapping }
        , _merge_branches{ merge_branches }
        , _algorithm{ algorithm }
        , _ghost_strategy{ strategy }
        , _kmer_size{ kmer_size }
        , _omega{ omega }
        , _filter{ filter }
        , _mu{ mu }
        , _num_threads{ num_threads }
        , _phylo_kmer_db{ kmer_size, omega, seq_type::name, i2l::io::to_newick(_original_tree)}
    {}

    void db_builder::run()
    {
        std::cout << "Construction parameters:" << std::endl <<
                  "\tSequence type: " << seq_type::name << std::endl <<
                  "\tk: " << _kmer_size << std::endl <<
                  "\tomega: " << _omega << std::endl <<
                  "\tKeep positions: " << (keep_positions ? "true" : "false") << std::endl << std::endl;

        /// Fill the tree index from the tree
        auto& index = _phylo_kmer_db.tree_index();
        index.reserve(_original_tree.get_node_count());
        for (const auto& node : visit_subtree(_original_tree.get_root()))
        {
            index.push_back(phylo_node::node_index{ node.get_num_nodes(), node.get_subtree_branch_length() });
        }

        /// The first stage of the algorithm: create a hashmap for every node group
        std::cout << "Building database: [stage 1 / 2]:" << std::endl;
        const auto& [group_ids, num_tuples, construction_time] = construct_group_hashmaps();
        std::cout << "Calculated " << num_tuples << " phylo-k-mers.\nCalculation time: " << construction_time
                  << "\n\n" << std::flush;

        /// The second stage of the algorithm: combine merge and filter
        std::cout << "Building database: [stage 2 / 2]:" << std::endl;
        const auto merge_time = merge(group_ids);
        //const auto merge_time = merge_filtered(group_ids);
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
    void normalize(std::vector<filter_value>& vec)
    {
        const auto min = std::min_element(vec.begin(), vec.end())->filter_score;
        const auto max = std::max_element(vec.begin(), vec.end())->filter_score;
        double sum = 0.0f;
        for (auto& elem : vec)
        {
            elem.filter_score = (elem.filter_score - min) / (max - min);
            sum += elem.filter_score;
        }

        for (auto& elem : vec)
        {
            elem.filter_score /= sum;
        }
    }

    void throw_if_positions()
    {
        #if defined(KEEP_POSITIONS)
            throw std::runtime_error("Positions are not supported in this version");
        #endif
    }

    unsigned long db_builder::merge(const std::vector<phylo_kmer::branch_type>& group_ids)
    {
        throw_if_positions();

        const auto begin = std::chrono::steady_clock::now();
        auto get_batch_db_name = [this](size_t batch_id) {
            return (fs::path(_working_directory) / fs::path{"hashmaps"} /
                    fs::path(std::to_string(batch_id) + ".ipk")).string();
        };

        std::cout << "Merge stage 1... ";
        /// Go over k-mer batches (ranges of k-mers)
        for (size_t batch_id = 0; batch_id < _num_batches; ++batch_id)
        {
            /// Merge all branch subdatabases for the current range of k-mers
            auto batch_db = ipk::merge_batch(_working_directory, group_ids, batch_id);

            const auto threshold = score_threshold(_omega, _kmer_size);
            auto filter = ipk::make_filter(_filter, _original_tree.get_node_count(),
                                           _working_directory, _num_batches, _mu, threshold);
            auto filter_values = filter->calc_filter_values(batch_db);
            normalize(filter_values);

            /// Sort filter values
            std::sort(filter_values.begin(), filter_values.end(),
                      [](const auto& a, const auto& b) { return a.filter_score > b.filter_score; });

            /// Save the order in which k-mer should be saved
            for (const auto& fv : filter_values)
            {
                batch_db.kmer_order.emplace_back(fv.key, (float)fv.filter_score);
            }

            /// Serialize the batch database
            i2l::save_uncompressed(batch_db, get_batch_db_name(batch_id));
        }
        std::cout << "OK" << std::endl;
        std::cout << "Merge stage 2...";

        /// Lazy loaders for every batch
        std::vector<batch_loader> batches;
        batches.reserve(_num_batches);
        // Priority queue for batch loaders
        std::priority_queue<batch_loader*, std::vector<batch_loader*>, batch_loader_compare> pq;
        for (size_t batch_id = 0; batch_id < _num_batches; ++batch_id)
        {
            batches.emplace_back(get_batch_db_name(batch_id));
            auto& loader = batches[batch_id];

            /// pre-load the first k-mer and push the loader to the queue
            if (loader.has_next())
            {
                loader.next();
                pq.push(&batches[batch_id]);
            }
        }

        // Lazy N-way merge of phylo-k-mers of all batches
        while (!pq.empty()) {
            batch_loader* loader = pq.top();
            pq.pop();

            if (auto& top = loader->current(); top.is_valid())
            {
                _phylo_kmer_db.insert_vector(top.key, std::move(top.entries));
            }

            // If the loader has more items, insert it back into the queue
            if (loader->has_next())
            {
                loader->next();
                pq.push(loader);
            }
        }
        std::cout << "OK" << std::endl;
        //_phylo_kmer_db.sort();

        const auto end = std::chrono::steady_clock::now();
        const auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

        std::cout << "Filtering and merge time: " << time << "\n\n" << std::flush;
        return time;
    }

    bool is_ghost(const phylo_node& node, ipk::ghost_strategy strategy)
    {
        const string& label = node.get_label();
        switch (strategy)
        {
            case ipk::ghost_strategy::INNER_ONLY:
                return boost::ends_with(label, "_X0");
            case ipk::ghost_strategy::OUTER_ONLY:
                return boost::ends_with(label, "_X1");
            default:
                return boost::ends_with(label, "_X0") || boost::ends_with(label, "_X1");
        }
    }

    /// \brief Returns a list of ghost node ids
    std::vector<std::string> get_ghost_ids(const phylo_tree& tree, ipk::ghost_strategy strategy)
    {
        std::vector<std::string> branch_ids;

        for (const auto& branch_node: tree)
        {
            if (is_ghost(branch_node, strategy))
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

        /// Filter and group ghost nodes
        const auto node_groups = group_ghost_ids(get_ghost_ids(_extended_tree, _ghost_strategy));

        /// Process branches in parallel. Results of the branch-and-bound algorithm are stored
        /// in a hash map for every group separately on disk.
        std::vector<phylo_kmer::branch_type> node_postorder_ids(node_groups.size());


        using namespace indicators;
        ProgressBar bar{
            option::BarWidth{60},
            option::Start{"["},
            option::Fill{"="},
            option::Lead{">"},
            option::Remainder{" "},
            option::End{"]"},
            option::PostfixText{"Computing phylo-k-mers"},
            option::ForegroundColor{Color::green},
            option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
            option::MaxProgress{node_groups.size()}
        };

        /*
        #pragma omp parallel for schedule(auto) reduction(+: count) num_threads(_num_threads) \
            shared(original_tree, node_postorder_ids, _matrix)
        */
        for (size_t i = 0; i < node_groups.size(); ++i)
        {
            const auto& node_group = node_groups[i];

            /// Having a label of a node in the extended tree, we need to find the corresponding node
            /// in the original tree. We take the first ghost node, because all of them correspond to
            /// the same original node
            const auto original_node_postorder_id = _extended_mapping.at(node_group[0]);
            node_postorder_ids[i] = original_node_postorder_id;


            /// Get sub-matrices of probabilities for the group
            const auto matrices = get_submatrices(node_group);

            /// Explore k-mers of the group and store results in a hash map
            const auto& [hash_maps, branch_count] = explore_group(matrices, original_node_postorder_id);

            /// Save the group hashmap on disk
            size_t index = 0;
            for (const auto& hash_map : hash_maps)
            {
                save_group_map(hash_map, get_group_map_file(_working_directory, original_node_postorder_id, index));
                ++index;
            }


            // update progress bar
            bar.set_option(option::PostfixText{std::to_string(i) + "/" + std::to_string(node_groups.size())});
            bar.tick();

            count += branch_count;
        }
        return { node_postorder_ids, count };
    }

    #ifdef KEEP_POSITIONS
    std::pair<std::vector<group_hash_map>, size_t> db_builder::explore_group(const proba_group& group, size_t postorder_id) const
    {
        (void)postorder_id;

        auto hash_maps = std::vector<group_hash_map>(_num_batches);
        size_t count = 0;

        const auto log_threshold = std::log10(i2l::score_threshold(_omega, _kmer_size));
        for (auto node_matrix_ref : group)
        {
            const auto& node_matrix = node_matrix_ref.get();
            for (const auto& window : to_windows(&node_matrix, _kmer_size))
            {
                auto alg = ipk::DCLA(window, _kmer_size);
                alg.run(log_threshold);
                for (const auto& kmer : alg.get_result())
                {
                    phylo_kmer positioned_kmer = { kmer.key, kmer.score,
                                                   static_cast<phylo_kmer::pos_type>(window.get_position()) };
                    put(hash_maps[kmer_batch(kmer.key, _num_batches)], positioned_kmer);
                    ++count;
                }
            }
        }

        return {std::move(hash_maps), count};
    }
    #else

    std::pair<std::vector<group_hash_map>, size_t> db_builder::explore_group(const proba_group& group, size_t postorder_id) const
    {
        (void)postorder_id;

        auto hash_maps = std::vector<group_hash_map>(_num_batches);

        size_t count = 0;
        const auto log_threshold = std::log10(score_threshold(_omega, _kmer_size));
        for (auto node_matrix_ref : group)
        {
            const auto& node_matrix = node_matrix_ref.get();
            for (const auto& window : to_windows(&node_matrix, _kmer_size))
            {
                auto alg = ipk::DCLA(window, _kmer_size);
                alg.run(log_threshold);

                for (const auto& kmer : alg.get_result())
                {
                    ipk::put(hash_maps[kmer_batch(kmer.key, _num_batches)], kmer);
                    ++count;
                }

            }
        }

        return { std::move(hash_maps), count };
    }

    #endif

};


namespace ipk
{
    phylo_kmer_db build(const string& working_directory,
                        const phylo_tree& original_tree, const phylo_tree& extended_tree,
                        const proba_matrix& matrix,
                        const ghost_mapping& mapping, const ar::mapping& ar_mapping, bool merge_branches,
                        ipk::algorithm algorithm, ipk::ghost_strategy strategy,
                        size_t kmer_size, i2l::phylo_kmer::score_type omega,
                        filter_type filter, double mu, size_t num_threads)
    {
        db_builder builder(working_directory,
                           original_tree, extended_tree,
                           matrix,
                           mapping, ar_mapping, merge_branches,
                           algorithm, strategy,
                           kmer_size, omega,
                           filter, mu, num_threads);
        builder.run();
        return std::move(builder._phylo_kmer_db);
    }
}
