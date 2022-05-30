#ifndef RAPPAS2_PK_COMPUTE_H
#define RAPPAS2_PK_COMPUTE_H


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
}


/*

#include <algorithm>
#include <cmath>
#include <xpas/phylo_kmer.h>
#include "common.h"
#include "matrix.h"

namespace xpas
{
    class dccw;

    class divide_and_conquer
    {
        friend class dccw;

    public:
        divide_and_conquer(const window& window, size_t k);
        void run(score_t omega);

        const map_t& get_map();

        const std::vector<unpositioned_phylo_kmer>& get_result() const;

        size_t get_num_kmers() const;

        void preprocess();

        std::vector<unpositioned_phylo_kmer> dc(score_t omega, size_t j, size_t h, score_t eps);
    private:


        score_t best_score(size_t j, size_t h);

        const window& _window;
        size_t _k;
        size_t _prefix_size;

        std::vector<unpositioned_phylo_kmer> _prefixes;
        std::vector<unpositioned_phylo_kmer> _suffixes;

        std::vector<score_t> _best_scores;

        std::vector<unpositioned_phylo_kmer> _result_list;
    };
}
 */
/*
class dccw
{
public:
    dccw(const window& window, std::vector<phylo_kmer>& prefixes, size_t k, score_t lookbehind, score_t lookahead);
    void run(score_t omega);

    const map_t& get_map();

    const std::vector<phylo_kmer>& get_result() const;

    size_t get_num_kmers() const;

    std::vector<phylo_kmer>&& get_suffixes();

    score_t get_best_suffix_score() const;

private:
    std::vector<phylo_kmer> dc(score_t omega, size_t j, size_t h, score_t eps);

    //void preprocess();

    //score_t best_score(size_t j, size_t h);

    const window& _window;
    size_t _k;
    size_t _prefix_size;

    // The second score bound for suffixes: the best suffix score of the next window
    score_t _lookahead;

    // The second score bound for prefixes: the best prefix score of the previous window
    score_t _lookbehind;

    std::vector<phylo_kmer>& _prefixes;
    std::vector<phylo_kmer> _suffixes;

    //std::vector<score_t> _best_scores;

    std::vector<phylo_kmer> _result_list;


    divide_and_conquer _dc;
};

*/


#endif //RAPPAS2_PK_COMPUTE_H
