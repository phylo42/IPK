#ifndef RAPPAS2_PK_COMPUTE_H
#define RAPPAS2_PK_COMPUTE_H

#include <vector>
#include "window.h"

namespace xpas
{
    enum class algorithm
    {
        // Branch-and-bound
        BB = 0,
        // Divide-and-conquer with no lookahead bound
        DC = 1,
        // Divide-and-conquer with the lookahead bound
        DCLA = 2,
        // Partition-based divide-and-conquer with Chained Windows
        DCCW = 3
    };


    class dccw;

    /// Divide-and-conquer with the lookahead trick
    class DCLA
    {
        friend class dccw;
    public:
        DCLA(const window& window, size_t k, phylo_kmer::score_type threshold);

        const std::vector<phylo_kmer>& get_result() const;

        std::vector<phylo_kmer> DC(size_t j, size_t h, phylo_kmer::score_type eps);
    private:

        void preprocess();

        phylo_kmer::score_type best_score(size_t j, size_t h);

        const window& _window;
        size_t _k;
        size_t _eps;
        size_t _prefix_size;

        std::vector<phylo_kmer> _prefixes;
        std::vector<phylo_kmer> _suffixes;

        std::vector<phylo_kmer::score_type> _best_scores;

        std::vector<phylo_kmer> _result_list;
    };


    class BB
    {
    public:
        BB(const window& window, size_t k, phylo_kmer::score_type eps);

        const std::vector<phylo_kmer>& get_result() const;
    private:
        void run(phylo_kmer::score_type eps);

        void bb(size_t i, size_t j, phylo_kmer::key_type prefix, phylo_kmer::score_type score, phylo_kmer::score_type eps);

        void preprocess();

        const window& _window;
        size_t _k;
        std::vector<phylo_kmer::score_type> _best_suffix_score;

        std::vector<phylo_kmer> _result_list;

        struct _mmer
        {
            phylo_kmer::key_type code;
            phylo_kmer::score_type score;
            unsigned short length;
        };
        std::vector<_mmer> _stack;
    };

}

#endif //RAPPAS2_PK_COMPUTE_H
