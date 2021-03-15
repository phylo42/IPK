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
    size_t num_entries;

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

kmer_filter::kmer_filter(std::string working_dir, size_t num_batches, double mu, phylo_kmer::score_type log_threshold,
                         xpas::score_model_type score_model)
    : _working_dir{ std::move(working_dir) }
    , _num_batches{ num_batches }
    , _mu{ mu }
    , _log_threshold{ log_threshold }
    , _score_model{ score_model }
{}

xpas::phylo_kmer::score_type logscore_to_score(xpas::phylo_kmer::score_type log_score)
{
    return std::min(std::pow(10, log_score), 1.0);
}

class batched_filter : public kmer_filter
{
public:
    batched_filter(size_t total_num_groups, std::string working_dir, size_t num_batches, double mu,
                   phylo_kmer::score_type log_threshold, xpas::score_model_type score_model)
        : kmer_filter{ std::move(working_dir), num_batches, mu, log_threshold, score_model },
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
            auto batch_db = merge_batch(_working_dir, group_ids, batch_idx, _score_model);
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

protected:
    size_t get_num_entries(const phylo_kmer_db& db) const
    {
        size_t num_entries = 0;
        for (const auto& [kmer, entries] : db)
        {
            num_entries += entries.size();
        }
        return num_entries;
    }
    virtual std::vector<filter_value> calc_filter_values(const phylo_kmer_db& db) const = 0;

    /// The total number of groups for which k-mers values are calculated
    size_t _total_num_groups;

    /// The set of k-mers taken by the filter
    std::unordered_set<phylo_kmer::key_type> _filtered;

};

class mif0_filter : public batched_filter
{
public:
    mif0_filter(size_t total_num_groups, std::string working_dir, size_t num_batches, double mu,
                phylo_kmer::score_type log_threshold, xpas::score_model_type score_model)
        : batched_filter{ total_num_groups, std::move(working_dir), num_batches, mu, log_threshold, score_model }
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
            const auto threshold = std::pow(10, _log_threshold);
            score_sum += static_cast<double>(_total_num_groups - entries.size()) * threshold;

            //std::cout << "\tTreshold: " << _threshold << std::endl;
            //std::cout << "\tS_w: " << score_sum << std::endl;

            /// s_wc / S_w, if s_wc == threshold
            const auto weighted_threshold = threshold / score_sum;
            const auto target_threshold = shannon(weighted_threshold);

            auto fv = filter_value{key, 0.0, entries.size()};
            auto HcBw1 = static_cast<double>(_total_num_groups) * target_threshold;
#ifdef KEEP_POSITIONS
            for (const auto& [branch, log_score, position] : entries)
#else
            for (const auto& [branch, log_score] : entries)
#endif
            {
                /// s_wc / S_w
                const auto weighted_score = logscore_to_score(log_score) / score_sum;
                const auto target_value = shannon(weighted_score);

                HcBw1 = HcBw1 - target_threshold + target_value;
                //std::cout << "FV: " << key << " " << xpas::decode_kmer(key, 5) << " " <<
                //                                                                      filter_stats[key] << std::endl;
            }

            //std::cout << "\tEntropy: " << filter_stats[key] << std::endl;
            const auto Hc = std::log2(_total_num_groups);

            /// MIs: Sw [ H(c) - H(c | B_w = 1) ] -> max
            ///   or equivalently
            ///      Sw [ H(c | B_w = 1) - H(c) ] -> min
            fv.filter_score = score_sum * (HcBw1 - Hc);
            //std::cout << "\t - MIF0: " << fv.filter_score << std::endl;

            filter_values.push_back(fv);
        }
        return filter_values;
    }
};


class mif1_filter : public batched_filter
{
public:
    mif1_filter(size_t total_num_groups, std::string working_dir, size_t num_batches, double mu,
                phylo_kmer::score_type log_threshold, xpas::score_model_type score_model)
        : batched_filter{ total_num_groups, std::move(working_dir), num_batches, mu, log_threshold, score_model }
    {}

    ~mif1_filter() noexcept override = default;

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
        const auto threshold = std::pow(10, _log_threshold);
        const auto N = _total_num_groups;

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
            Sw += static_cast<double>(N - entries.size()) * threshold;

            //std::cout << "\tSW = " << Sw << std::endl;

            /// Calculate the Mutual Information
            /// MI = H(C) - H(C|Bw) =
            ///    = H(C) - [P(Bw=1)H(C|Bw=1) + P(Bw=0)H(C|Bw=0)] =
            ///    = A - (B + C)
            const auto A = std::log2(_total_num_groups);

            /// P(Bw == 1) = s_wc / Sw, if s_wc == threshold
            const auto PBw1_threshold = threshold / Sw;
            /// P(Bw == 0) = (1 - s_wc) / (N - Sw), if s_wc == threshold
            const auto PBw0_threshold = (1 - threshold) / (N - Sw);

            /// Calculate B and C: starting with the missing scores
            double HcBw1 = static_cast<double>(N) * shannon(PBw1_threshold);
            double HcBw0 = static_cast<double>(N) * shannon(PBw0_threshold);
            auto fv = filter_value{key, 0.0, entries.size()};
#ifdef KEEP_POSITIONS
            for (const auto& [branch, log_score, position] : entries)
#else
            for (const auto& [branch, log_score] : entries)
#endif
            {
                const auto swc = logscore_to_score(log_score);

                /// P(c|Bw == 1) = s_wc / Sw
                const auto PcBw1 = swc / Sw;
                /// P(c|Bw == 0) = (1 - s_wc) / (N - Sw)
                const auto PcBw0 = (1 - swc) / (N - Sw);

                HcBw1 = HcBw1 - shannon(PBw1_threshold) + shannon(PcBw1);
                HcBw0 = HcBw0 - shannon(PBw0_threshold) + shannon(PcBw0);
            }
            /// P(Bw=1)
            const auto PBw1 = Sw / N;
            const auto B = PBw1 * HcBw1;
            /// P(Bw=0)
            const auto PBw0 = (N - Sw) / N;
            const auto C = PBw0 * HcBw0;


            /// A - (B + C) -> max
            /// which is A - B - C -> max
            /// which is -A + B + C -> min
            fv.filter_score = - A + B + C;
            filter_values.push_back(fv);
            //std::cout << "\tA: " << A << std::endl;
            //std::cout << "\tB: " << B << std::endl;
            //std::cout << "\tC: " << C << std::endl;
            //std::cout << "\tMIF1: " << fv.filter_score << std::endl;
        }

        return filter_values;
    }
};

class random_filter : public batched_filter
{
public:
    random_filter(size_t total_num_groups, std::string working_dir, size_t num_batches, double mu,
                  phylo_kmer::score_type log_threshold, xpas::score_model_type score_model)
        : batched_filter{ total_num_groups, std::move(working_dir), num_batches, mu, log_threshold, score_model }
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

class no_filter : public kmer_filter
{
public:
    no_filter(std::string working_dir, size_t num_batches, double mu,
              phylo_kmer::score_type log_threshold, xpas::score_model_type score_model)
        : kmer_filter{ std::move(working_dir), num_batches, mu, log_threshold, score_model }
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
                                               double mu,
                                               phylo_kmer::score_type log_threshold, xpas::score_model_type score_model)
{
    if (filter == filter_type::no_filter || mu == 1.0)
    {
        return std::make_unique<no_filter>(std::move(working_dir), num_batches, mu,
                                           log_threshold, score_model);
    }
    else if (filter == filter_type::mif0)
    {
        return std::make_unique<mif0_filter>(total_num_nodes, std::move(working_dir), num_batches, mu,
                                             log_threshold, score_model);
    }
    else if (filter == filter_type::mif1)
    {
        return std::make_unique<mif1_filter>(total_num_nodes, std::move(working_dir), num_batches, mu,
                                             log_threshold, score_model);
    }
    else if (filter == filter_type::random)
    {
        return std::make_unique<random_filter>(total_num_nodes, std::move(working_dir), num_batches, mu,
                                               log_threshold, score_model);
    }
    else
    {
        throw std::runtime_error("Error: Unsupported filter type.");
    }
}