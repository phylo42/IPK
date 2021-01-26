#include "filter.h"
#include "branch_group.h"
#include <xpas/phylo_kmer_db.h>
#include <boost/filesystem.hpp>
#include <vector>
#include <queue>
#include <unordered_set>
#include <iostream>
#include <random>

using namespace xpas;
namespace fs = boost::filesystem;


/// A pair storing the information needed to evaluate a phylo k-mer by
/// the filtering function. "Value" is the value calculated by the filter
/// which is used to compare the informativeness of the k-mer against the others
struct filter_value
{
    phylo_kmer::key_type key;
    double filter_score;

    bool operator<(const filter_value& rhs) const
    {
        return filter_score < rhs.filter_score;
    }

    bool operator>(const filter_value& rhs) const
    {
        return filter_score > rhs.filter_score;
    }
};

std::string xpas::get_fvs_file(const std::string& working_dir, size_t batch_idx)
{
    return (fs::path{xpas::get_groups_dir(working_dir)} / fs::path(std::to_string(batch_idx) + ".fvs")).string();
}

kmer_filter::kmer_filter(std::string working_dir, size_t num_batches, double mu, phylo_kmer::score_type threshold)
    : _working_dir{ std::move(working_dir) }, _num_batches{ num_batches }, _mu{ mu }, _threshold{ threshold }
{}

xpas::phylo_kmer::score_type logscore_to_score(xpas::phylo_kmer::score_type log_score)
{
    return std::min(std::pow(10, log_score), 1.0);
}

class entropy_filter : public kmer_filter
{
public:
    entropy_filter(size_t total_num_groups, std::string working_dir, size_t num_batches, double mu, phylo_kmer::score_type threshold)
        : kmer_filter{ std::move(working_dir), num_batches, mu, threshold },
        _total_num_groups{ total_num_groups }
    {}

    ~entropy_filter() noexcept override = default;

    void filter(const std::vector<phylo_kmer::branch_type>& group_ids) override
    {
        std::cout << "Filtering (minimal conditional entropy)..." << std::endl;

        /// Create Filter Value files for every batch
        std::vector<std::vector<filter_value>> batch_filter_values;

        size_t total_num_kmers = 0;
        for (size_t batch_idx = 0; batch_idx < _num_batches; ++batch_idx)
        {
            /// Merge hashmaps of the same batch
            auto batch_db = merge_batch(_working_dir, group_ids, batch_idx);
            total_num_kmers += batch_db.size();

            /// Calculate filter values for the batch
            auto filter_values = calc_filter_values(batch_db);

            //FIXME: do the partial sort to the (mu * total_num_kmers)th element
            std::sort(filter_values.begin(), filter_values.end());
            batch_filter_values.push_back(std::move(filter_values));
        }

        /// The vector of indexes of top kmers for every batch,
        /// needed for the merge
        std::vector<size_t> batch_idx(_num_batches, 0);

        /// Initialize a min-heap for the merge algorithm
        std::priority_queue<filter_value, std::vector<filter_value>, std::greater<>> heap;
        for (const auto& batch_fv : batch_filter_values)
        {
            heap.push(batch_fv[0]);
        }

        size_t final_num_kmers = _mu * total_num_kmers;
        for (size_t i = 0; i < final_num_kmers; ++i)
        {
            /// Mark the top k-mer as good
            auto best_fv = heap.top();
            heap.pop();
            _filtered.insert(best_fv.key);
            std::cout << "TAKE " << best_fv.key << " " << xpas::decode_kmer(best_fv.key, 5) << " " <<
                best_fv.filter_score << std::endl;

            /// take the next k-mer value from this batch
            size_t best_kmer_index = kmer_batch(best_fv.key, _num_batches);
            const auto& best_kmer_batch = batch_filter_values[best_kmer_index];
            batch_idx[best_kmer_index]++;
            if (batch_idx[best_kmer_index] < best_kmer_batch.size())
            {
                /// and push it to the heap
                heap.push(best_kmer_batch[batch_idx[best_kmer_index]]);
            }
        }
    }

    bool is_good(phylo_kmer::key_type key) const noexcept override
    {
        return _filtered.find(key) != _filtered.end();
    }

private:
    static double shannon(double x)
    {
        return - x * std::log2(x);
    }

    std::vector<filter_value> calc_filter_values(const phylo_kmer_db& db) const
    {
        std::vector<filter_value> filter_values;

        hash_map<phylo_kmer::key_type, double> filter_stats;
        for (const auto& [key, entries] : db)
        {
            /// calculate the score sum to normalize scores
            double score_sum = 0;
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
            }
#endif

            /// do not forget the branches that are not stored in the database,
            /// they suppose to have the threshold score
            score_sum += static_cast<double>(_total_num_groups - entries.size()) * _threshold;

            /// Entropy
            const auto weighted_threshold = _threshold / score_sum;
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
                    filter_stats[key] = static_cast<double>(_total_num_groups) * target_threshold;
                }
                filter_stats[key] = filter_stats[key] - target_threshold + target_value;
                //std::cout << "FV: " << key << " " << xpas::decode_kmer(key, 5) << " " <<
                //                                                                      filter_stats[key] << std::endl;
            }
        }
        /// copy entropy values into a vector
        filter_values.reserve(filter_stats.size());
        for (const auto& [key, filter_score] : filter_stats)
        {
            filter_values.push_back({ key, filter_score });
        }
        return filter_values;
    }

    /// The total number of groups for which k-mers values are calculated
    size_t _total_num_groups;

    /// The set of k-mers taken by the filter
    std::unordered_set<phylo_kmer::key_type> _filtered;
};


class random_filter : public kmer_filter
{
public:
    random_filter(std::string working_dir, size_t num_batches, double mu, phylo_kmer::score_type threshold)
        : kmer_filter{ std::move(working_dir), num_batches, mu, threshold }
    {}

    ~random_filter() noexcept = default;

    void filter(const std::vector<phylo_kmer::branch_type>& group_ids)
    {
        std::default_random_engine generator(42);
        std::uniform_real_distribution<double> distribution(0, 1);

        size_t total_num_kmers = 0;
        for (size_t batch_idx = 0; batch_idx < _num_batches; ++batch_idx)
        {
            /// Merge hashmaps of ghost nodes of the same batch
            auto batch_db = merge_batch(_working_dir, group_ids, batch_idx);
            total_num_kmers += batch_db.size();

            /// update the map randomly
            for (const auto& [key, entries] : batch_db)
            {
                if (distribution(generator) <= _mu)
                {
                    _filtered.insert(key);
                }
            }
        }
    }

    bool is_good(phylo_kmer::key_type key) const noexcept override
    {
        return _filtered.find(key) != _filtered.end();
    }
private:

    /// The set of k-mers taken by the filter
    std::unordered_set<phylo_kmer::key_type> _filtered;
};

class no_filter : public kmer_filter
{
public:
    no_filter(std::string working_dir, size_t num_batches, double mu, phylo_kmer::score_type threshold)
        : kmer_filter{ std::move(working_dir), num_batches, mu, threshold }
    {}

    ~no_filter() noexcept = default;

    void filter(const std::vector<phylo_kmer::branch_type>& group_ids)
    {
        (void)group_ids;
    }

    bool is_good(phylo_kmer::key_type key) const noexcept override
    {
        (void)key;
        return true;
    }

};


std::unique_ptr<kmer_filter> xpas::make_filter(xpas::filter_type filter,
                                               size_t total_num_nodes,
                                               std::string working_dir, size_t num_batches,
                                               double mu, phylo_kmer::score_type threshold)
{
    if (filter == filter_type::entropy)
    {
        return std::make_unique<entropy_filter>(total_num_nodes, std::move(working_dir), num_batches, mu, threshold);
    }
    else if (filter == filter_type::random)
    {
        return std::make_unique<random_filter>(std::move(working_dir), num_batches, mu, threshold);
    }
    else if (filter == filter_type::no_filter)
    {
        return std::make_unique<no_filter>(std::move(working_dir), num_batches, mu, threshold);
    }
    else
    {
        throw std::runtime_error("Error: Unsupported filter type.");
    }
}