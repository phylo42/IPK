#include "filter.h"
#include "branch_group.h"
#include <xpas/phylo_kmer_db.h>
#include <xpas/version.h>
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
        for (size_t batch_idx = 0; batch_idx < _num_batches; ++batch_idx)
        {
            /// Merge hashmaps of the same batch
            auto batch_db = merge_batch(_working_dir, group_ids, batch_idx);
            total_num_kmers += batch_db.size();

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

        size_t final_num_kmers = _mu * total_num_kmers;
        for (size_t i = 0; i < final_num_kmers; ++i)
        {
            /// Mark the top k-mer as good
            auto best_fv = heap.top();
            heap.pop();
            _filtered.insert(best_fv.key);

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

protected:
    virtual std::vector<filter_value> calc_filter_values(const phylo_kmer_db& db) const = 0;

    /// The total number of groups for which k-mers values are calculated
    size_t _total_num_groups;

    /// The set of k-mers taken by the filter
    std::unordered_set<phylo_kmer::key_type> _filtered;

};

class entropy_filter : public batched_filter
{
public:
    entropy_filter(size_t total_num_groups, std::string working_dir, size_t num_batches, double mu,
                   phylo_kmer::score_type threshold)
        : batched_filter{ total_num_groups, std::move(working_dir), num_batches, mu, threshold }
    {}

    ~entropy_filter() noexcept override = default;

    void filter(const std::vector<phylo_kmer::branch_type>& group_ids) override
    {
        std::cout << "Filtering (minimal conditional entropy)..." << std::endl;
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
            }
#else
            for (const auto& [branch, log_score] : entries)
            {
                (void)branch;
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

};

class mis_filter : public batched_filter
{
public:
    mis_filter(size_t total_num_groups, std::string working_dir, size_t num_batches, double mu,
                   phylo_kmer::score_type threshold)
        : batched_filter{ total_num_groups, std::move(working_dir), num_batches, mu, threshold }
    {}

    ~mis_filter() noexcept override = default;

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

        hash_map<phylo_kmer::key_type, double> filter_stats;
        for (const auto& [key, entries] : db)
        {
            /// calculate the score sum to normalize scores7
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

                //std::cout << "\t" << logscore_to_score(log_score) << "," << std::endl;
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

#ifdef KEEP_POSITIONS
            for (const auto& [branch, log_score, position] : entries)
#else
            for (const auto& [branch, log_score] : entries)
#endif
            {
                /// s_wc / S_w
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

            //std::cout << "\tEntropy: " << filter_stats[key] << std::endl;

            /// MIs = S_w / N * H(c | B_w = 1)
            filter_stats[key] = score_sum / _total_num_groups * filter_stats[key];

            //std::cout << "\tMIS: " << filter_stats[key] << std::endl;
        }
        /// copy entropy values into a vector
        filter_values.reserve(filter_stats.size());
        for (const auto& [key, filter_score] : filter_stats)
        {
            filter_values.push_back({ key, filter_score });
        }
        return filter_values;
    }
};


class mif_filter : public batched_filter
{
public:
    mif_filter(size_t total_num_groups, std::string working_dir, size_t num_batches, double mu,
               phylo_kmer::score_type threshold)
        : batched_filter{ total_num_groups, std::move(working_dir), num_batches, mu, threshold }
    {}

    ~mif_filter() noexcept override = default;

    void filter(const std::vector<phylo_kmer::branch_type>& group_ids) override
    {
        std::cout << "Filtering (maximal mutual information, full version)..." << std::endl;
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
        const auto N = _total_num_groups;

        hash_map<phylo_kmer::key_type, double> filter_stats;
        for (const auto& [key, entries] : db)
        {
            //std::cout << key << std::endl;

            /// calculate the score sum to normalize scores
            double Sw = 0;
#ifdef KEEP_POSITIONS
            for (const auto& [branch, log_score, position] : entries)
            {
                (void)branch;
                (void)position;
                Sw += logscore_to_score(log_score);
            }
#else
            for (const auto& [branch, log_score] : entries)
            {
                (void)branch;
                Sw += logscore_to_score(log_score);

                //std::cout << "\t" << logscore_to_score(log_score) << "," << std::endl;
            }
#endif

            /// do not forget the branches that are not stored in the database,
            /// they suppose to have the threshold score
            Sw += static_cast<double>(N - entries.size()) * _threshold;

            //std::cout << "\tSW = " << Sw << std::endl;

            /// Calculate the Mutual Information
            /// MI = H(C) - H(C|Bw) =
            ///    = H(C) - [P(Bw=1)H(C|Bw=1) + P(Bw=0)H(C|Bw=0)] =
            ///    = A - (B + C)

            /// P(Bw == 1) = s_wc / Sw, if s_wc == threshold
            const auto PBw1_threshold = _threshold / Sw;
            /// P(Bw == 0) = (1 - s_wc) / (N - Sw), if s_wc == threshold
            const auto PBw0_threshold = (1 - _threshold) / (N - Sw);

            /// Calculate B and C: starting with the missing scores
            double HcBw1 = static_cast<double>(N) * shannon(PBw1_threshold);
            double HcBw0 = static_cast<double>(N) * shannon(PBw0_threshold);

#ifdef KEEP_POSITIONS
            for (const auto& [branch, log_score, position] : entries)
#else
            for (const auto& [branch, log_score] : entries)
#endif
            {
                const auto swc = logscore_to_score(log_score);

                /// P(Bw == 1) = s_wc / Sw
                const auto PBw1 = swc / Sw;
                /// P(Bw == 0) = (1 - s_wc) / (N - Sw)
                const auto PBw0 = (1 - swc) / (N - Sw);

                HcBw1 = HcBw1 - shannon(PBw1_threshold) + shannon(PBw1);
                HcBw0 = HcBw0 - shannon(PBw0_threshold) + shannon(PBw0);
            }

            /// B = P(Bw=1) * H(c|Bw=1)
            const auto b = Sw / N * HcBw1;

            /// C = P(Bw=0) * H(c|Bw=0)
            const auto c = (N - Sw) / N * HcBw0;

            /// Since A is const for all k-mers,
            /// maximizing A - (B + C) is the same as minimizing B + C
            filter_stats[key] = b + c;

            //std::cout << "\tB: " << b << std::endl;
            //std::cout << "\tC: " << c << std::endl;
            //std::cout << "\tMIF: " << filter_stats[key] << std::endl;
        }

        /// copy entropy values into a vector
        filter_values.reserve(filter_stats.size());
        for (const auto& [key, filter_score] : filter_stats)
        {
            filter_values.push_back({ key, filter_score });
        }
        return filter_values;
    }
};


class random_filter : public kmer_filter
{
public:
    random_filter(std::string working_dir, size_t num_batches, double mu, phylo_kmer::score_type threshold)
        : kmer_filter{ std::move(working_dir), num_batches, mu, threshold }
    {}

    ~random_filter() noexcept override = default;

    void filter(const std::vector<phylo_kmer::branch_type>& group_ids) override
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

    ~no_filter() noexcept override = default;

    void filter(const std::vector<phylo_kmer::branch_type>& group_ids) override
    {
        (void)group_ids;
    }

    [[nodiscard]]
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
    if (filter == filter_type::no_filter || mu == 1.0)
    {
        return std::make_unique<no_filter>(std::move(working_dir), num_batches, mu, threshold);
    }
    else if (filter == filter_type::entropy)
    {
        return std::make_unique<entropy_filter>(total_num_nodes, std::move(working_dir), num_batches, mu, threshold);
    }
    else if (filter == filter_type::mis)
    {
        return std::make_unique<mis_filter>(total_num_nodes, std::move(working_dir), num_batches, mu, threshold);
    }
    else if (filter == filter_type::mif)
    {
        return std::make_unique<mif_filter>(total_num_nodes, std::move(working_dir), num_batches, mu, threshold);
    }
    else if (filter == filter_type::random)
    {
        return std::make_unique<random_filter>(std::move(working_dir), num_batches, mu, threshold);
    }
    else
    {
        throw std::runtime_error("Error: Unsupported filter type.");
    }
}