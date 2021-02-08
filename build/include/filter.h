#ifndef XPAS_FILTER_H
#define XPAS_FILTER_H

#include <xpas/phylo_kmer.h>
#include <memory>

namespace xpas
{
    enum class filter_type
    {
        no_filter,
        entropy,
        random,
        mis,
        mif
    };

    /// Get the .fvs file of the k-mer batch (filter value stats)
    std::string get_fvs_file(const std::string& working_dir, size_t batch_idx);

    /// Determines if the given k-mer has to be filtered out or not.
    class kmer_filter
    {
    public:
        kmer_filter(std::string working_dir, size_t num_batches, double mu, phylo_kmer::score_type threshold);

        virtual ~kmer_filter() noexcept = default;

        virtual void filter(const std::vector<phylo_kmer::branch_type>& group_ids) = 0;

        [[nodiscard]]
        virtual bool is_good(phylo_kmer::key_type key) const noexcept = 0;
    protected:
        std::string _working_dir;
        size_t _num_batches;
        double _mu;
        phylo_kmer::score_type _threshold;
    };

    std::unique_ptr<kmer_filter> make_filter(xpas::filter_type filter,
                                             size_t total_num_nodes,
                                             std::string working_dir, size_t num_batches,
                                             double mu, phylo_kmer::score_type threshold);

}


#endif //XPAS_FILTER_H
