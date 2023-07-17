#include "filter.h"
#include "branch_group.h"
#include <i2l/phylo_kmer_db.h>
#include <i2l/version.h>
#include <boost/filesystem.hpp>
#include <vector>
#include <queue>
#include <unordered_set>
#include <iostream>
#include <random>

using namespace i2l;
using namespace ipk;
namespace fs = boost::filesystem;


bool ipk::filter_value::operator<(const ipk::filter_value& rhs) const
{
    return filter_score < rhs.filter_score;
}

bool ipk::filter_value::operator>(const ipk::filter_value& rhs) const
{
    return filter_score > rhs.filter_score;
}

kmer_filter::kmer_filter(std::string working_dir, size_t num_batches, double mu, phylo_kmer::score_type threshold)
    : _working_dir{ std::move(working_dir) }, _num_batches{ num_batches }, _mu{ mu }, _threshold{ threshold }
{}

i2l::phylo_kmer::score_type logscore_to_score(i2l::phylo_kmer::score_type log_score)
{
    return std::min(std::pow(10, log_score), 1.0);
}

class batched_filter : public kmer_filter
{
public:
    batched_filter(size_t total_num_groups, std::string working_dir, size_t num_batches, double mu, phylo_kmer::score_type threshold)
        : kmer_filter{ std::move(working_dir), num_batches, mu, threshold },
        _total_num_groups{ total_num_groups }
    {}

    ~batched_filter() noexcept override = default;

    void filter(const std::vector<phylo_kmer::branch_type>& group_ids) override
    {
        /// Create Filter Value files for every batch
        std::vector<std::vector<filter_value>> batch_filter_values;

        size_t total_num_kmers = 0;

        /// num of pairs { branch -> score }
        size_t total_num_entries = 0;

        for (size_t batch_idx = 0; batch_idx < _num_batches; ++batch_idx)
        {
            /// Merge hashmaps of the same batch
            auto batch_db = merge_batch(_working_dir, group_ids, batch_idx);
            total_num_kmers += batch_db.size();
            total_num_entries += get_num_entries(batch_db);

            /// Calculate filter values for the batch
            auto filter_values = calc_filter_values(batch_db);

            /// FIXME: do the partial sort to the (mu * total_num_kmers)th element
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

        size_t max_entries_allowed = _mu * total_num_entries;
        size_t num_entries = 0;
        while (num_entries < max_entries_allowed)
        {
            /// Mark the top k-mer as good
            auto best_fv = heap.top();
            heap.pop();
            _filtered.insert(best_fv.key);
            num_entries += best_fv.num_entries;

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

    virtual std::vector<filter_value> calc_filter_values(const phylo_kmer_db& db) const = 0;

protected:
    size_t get_num_entries(const phylo_kmer_db& db) const
    {
        size_t num_entries = 0;
        for (const auto& [kmer, entries] : db)
        {
            (void)kmer;
            num_entries += entries.size();
        }
        return num_entries;
    }

    /// The total number of groups for which k-mers values are calculated
    size_t _total_num_groups;

    /// The set of k-mers taken by the filter
    std::unordered_set<phylo_kmer::key_type> _filtered;

};


class mif0_filter : public batched_filter
{
public:
    mif0_filter(size_t total_num_groups, std::string working_dir, size_t num_batches, double mu,
                phylo_kmer::score_type threshold)
        : batched_filter{ total_num_groups, std::move(working_dir), num_batches, mu, threshold }
    {}

    ~mif0_filter() noexcept override = default;

    void filter(const std::vector<phylo_kmer::branch_type>& group_ids) override
    {
        std::cout << "Filtering (maximal mutual information, short version)..." << std::endl;
        batched_filter::filter(group_ids);
    }

private:
    static double shannon(double x)
    {
        return - x * std::log2(x);
    }

    std::vector<filter_value> calc_filter_values(const phylo_kmer_db& db) const override
    {
        std::vector<filter_value> filter_values;
        for (const auto& [key, entries] : db)
        {
            /// calculate the score sum to normalize scores
            /// i.e. S_w by the notation
            double score_sum = 0;
#ifdef KEEP_POSITIONS
            for (const auto& [branch, log_score, position] : entries)
            {
                (void)branch;
                (void)position;
                score_sum += logscore_to_score(log_score);
            }
#else
            //std::cout << key << std::endl;
            for (const auto& [branch, log_score] : entries)
            {
                (void)branch;
                score_sum += logscore_to_score(log_score);
                //std::cout << "\t" << logscore_to_score(log_score) << ",";
            }
#endif

            /// do not forget the branches that are not stored in the database,
            /// they suppose to have the threshold score
            score_sum += static_cast<double>(_total_num_groups - entries.size()) * _threshold;

            //std::cout << "\tTreshold: " << _threshold << std::endl;
            //std::cout << "\tS_w: " << score_sum << std::endl;

            /// s_wc / S_w, if s_wc == threshold
            const auto weighted_threshold = _threshold / score_sum;
            const auto target_threshold = shannon(weighted_threshold);

            auto fv = filter_value{key, 0.0, entries.size()};
            auto HcBw1 = static_cast<double>(_total_num_groups) * target_threshold;
#ifdef KEEP_POSITIONS
            for (const auto& [branch, log_score, position] : entries)
            {
                (void)position;
#else
            for (const auto& [branch, log_score] : entries)
            {
#endif
                (void)branch;
                /// s_wc / S_w
                const auto weighted_score = logscore_to_score(log_score) / score_sum;
                const auto target_value = shannon(weighted_score);

                HcBw1 = HcBw1 - target_threshold + target_value;
            }

            //std::cout << "\tEntropy: " << filter_stats[key] << std::endl;
            const auto Hc = std::log2(_total_num_groups);

            /// MIs: Sw [ H(c) - H(c | B_w = 1) ] -> max
            ///   or equivalently
            ///      Sw [ H(c | B_w = 1) - H(c) ] -> min
            fv.filter_score = score_sum * (HcBw1 - Hc);
            //std::cout << fv.filter_score << ",";

            filter_values.push_back(fv);
        }
        return filter_values;
    }
};

class random_filter : public batched_filter
{
public:
    random_filter(size_t total_num_groups, std::string working_dir,
                  size_t num_batches, double mu, phylo_kmer::score_type threshold)
        : batched_filter{ total_num_groups, std::move(working_dir), num_batches, mu, threshold }
    {}

    ~random_filter() noexcept override = default;

    void filter(const std::vector<phylo_kmer::branch_type>& group_ids) override
    {
        std::cout << "Filtering (random)..." << std::endl;
        batched_filter::filter(group_ids);
    }

private:
    std::vector<filter_value> calc_filter_values(const phylo_kmer_db& db) const override
    {
        std::vector<filter_value> filter_values;

        std::default_random_engine generator(42);
        std::uniform_real_distribution<double> distribution(0, 1);

        for (const auto& [key, entries] : db)
        {
            auto fv = filter_value{ key, distribution(generator), entries.size() };
            filter_values.push_back(fv);
        }
        return filter_values;
    }
};

std::unique_ptr<kmer_filter> ipk::make_filter(ipk::filter_type filter,
                                              size_t total_num_nodes,
                                              std::string working_dir, size_t num_batches,
                                              double mu, phylo_kmer::score_type threshold)
{
    if (filter == filter_type::mif0)
    {
        return std::make_unique<mif0_filter>(total_num_nodes, std::move(working_dir), num_batches, mu, threshold);
    }
    else if (filter == filter_type::random)
    {
        return std::make_unique<random_filter>(total_num_nodes, std::move(working_dir), num_batches, mu, threshold);
    }
    else
    {
        throw std::runtime_error("Error: Unsupported filter type.");
    }
}