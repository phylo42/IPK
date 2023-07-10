#ifndef XPAS_FILTER_H
#define XPAS_FILTER_H

#include <i2l/phylo_kmer.h>
#include <i2l/phylo_kmer_db.h>
#include <memory>

namespace ipk
{
    using i2l::phylo_kmer;

    enum class filter_type
    {
        random,
        mif0,
    };

    /// A pair storing the information needed to evaluate a phylo k-mer by
    /// the filtering function. "Value" is the value calculated by the filter
    /// which is used to compare the informativeness of the k-mer against the others
    struct filter_value
    {
        phylo_kmer::key_type key;
        double filter_score;
        size_t num_entries;

        bool operator<(const filter_value& rhs) const;
        bool operator>(const filter_value& rhs) const;
    };

        /// Determines if the given k-mer has to be filtered out or not.
    class kmer_filter
    {
    public:
        kmer_filter(std::string working_dir, size_t num_batches, double mu, phylo_kmer::score_type threshold);

        virtual ~kmer_filter() noexcept = default;

        virtual void filter(const std::vector<phylo_kmer::branch_type>& group_ids) = 0;

        [[nodiscard]]
        virtual std::vector<filter_value> calc_filter_values(const i2l::phylo_kmer_db& db) const = 0;

        [[nodiscard]]
        virtual bool is_good(phylo_kmer::key_type key) const noexcept = 0;
    protected:
        std::string _working_dir;
        size_t _num_batches;
        double _mu;
        phylo_kmer::score_type _threshold;
    };

    std::unique_ptr<kmer_filter> make_filter(ipk::filter_type filter,
                                             size_t total_num_nodes,
                                             std::string working_dir, size_t num_batches,
                                             double mu, phylo_kmer::score_type threshold);

}


#endif //XPAS_FILTER_H
