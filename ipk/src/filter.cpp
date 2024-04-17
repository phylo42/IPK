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

kmer_filter::kmer_filter(std::string working_dir, size_t num_batches, phylo_kmer::score_type threshold)
    : _working_dir{ std::move(working_dir) }, _num_batches{ num_batches }, _threshold{ threshold }
{}

i2l::phylo_kmer::score_type logscore_to_score(i2l::phylo_kmer::score_type log_score)
{
    return std::min(std::pow(10, log_score), 1.0);
}

class batched_filter : public kmer_filter
{
public:
    batched_filter(size_t total_num_groups, std::string working_dir,
                   size_t num_batches, phylo_kmer::score_type threshold)
        : kmer_filter{ std::move(working_dir), num_batches, threshold },
        _total_num_groups{ total_num_groups }
    {}

    ~batched_filter() noexcept override = default;

    virtual std::vector<kmer_fv> calc_filter_values(const phylo_kmer_db& db) const = 0;

protected:
    /// The total number of groups for which k-mers values are calculated
    size_t _total_num_groups;
};


class mif0_filter : public batched_filter
{
public:
    mif0_filter(size_t total_num_groups, std::string working_dir, size_t num_batches,
                phylo_kmer::score_type threshold)
        : batched_filter{ total_num_groups, std::move(working_dir), num_batches, threshold }
    {}

    ~mif0_filter() noexcept override = default;

private:
    static double shannon(double x)
    {
        return - x * std::log2(x);
    }

    std::vector<kmer_fv> calc_filter_values(const phylo_kmer_db& db) const override
    {
        std::vector<kmer_fv> filter_values;
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
            }
#endif

            /// do not forget the branches that are not stored in the database,
            /// they suppose to have the threshold score
            score_sum += static_cast<double>(_total_num_groups - entries.size()) * _threshold;

            /// s_wc / S_w, if s_wc == threshold
            const auto weighted_threshold = _threshold / score_sum;
            const auto target_threshold = shannon(weighted_threshold);

            auto fv = kmer_fv{ key, 0.0 };
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

            const auto Hc = std::log2(_total_num_groups);

            /// MIs: Sw [ H(c) - H(c | B_w = 1) ] -> max
            ///   or equivalently
            ///      Sw [ H(c | B_w = 1) - H(c) ] -> min
            fv.filter_value = score_sum * (HcBw1 - Hc);
            filter_values.push_back(fv);
        }
        return filter_values;
    }
};

class random_filter : public batched_filter
{
public:
    random_filter(size_t total_num_groups, std::string working_dir,
                  size_t num_batches, phylo_kmer::score_type threshold)
        : batched_filter{ total_num_groups, std::move(working_dir), num_batches, threshold }
    {}

    ~random_filter() noexcept override = default;

private:
    std::vector<kmer_fv> calc_filter_values(const phylo_kmer_db& db) const override
    {
        std::vector<kmer_fv> filter_values;

        std::default_random_engine generator(42);
        std::uniform_real_distribution<double> distribution(0, 1);

        for (const auto& [key, entries] : db)
        {
            auto fv = kmer_fv{ key, (float)distribution(generator) };
            filter_values.push_back(fv);
        }
        return filter_values;
    }
};

std::unique_ptr<kmer_filter> ipk::make_filter(ipk::filter_type filter,
                                              size_t total_num_nodes,
                                              std::string working_dir, size_t num_batches,
                                              phylo_kmer::score_type threshold)
{
    if (filter == filter_type::mif0)
    {
        return std::make_unique<mif0_filter>(total_num_nodes, std::move(working_dir), num_batches, threshold);
    }
    else if (filter == filter_type::random)
    {
        return std::make_unique<random_filter>(total_num_nodes, std::move(working_dir), num_batches, threshold);
    }
    else
    {
        throw std::runtime_error("Error: Unsupported filter type.");
    }
}