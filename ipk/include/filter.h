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

    /// Determines if the given k-mer has to be filtered out or not.
    class kmer_filter
    {
    public:
        kmer_filter(std::string working_dir, size_t num_batches, phylo_kmer::score_type threshold);

        virtual ~kmer_filter() noexcept = default;

        [[nodiscard]]
        virtual std::vector<i2l::kmer_fv> calc_filter_values(const i2l::phylo_kmer_db& db) const = 0;

    protected:
        std::string _working_dir;
        size_t _num_batches;
        phylo_kmer::score_type _threshold;
    };

    std::unique_ptr<kmer_filter> make_filter(ipk::filter_type filter,
                                             size_t total_num_nodes,
                                             std::string working_dir, size_t num_batches,
                                             phylo_kmer::score_type threshold);

}


#endif //XPAS_FILTER_H
