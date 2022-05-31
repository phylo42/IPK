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


    namespace impl
    {
        /* Divide-and-conquer with the lookahead technique
        class DCLA_iterator
        {
        public:
            DCLA_iterator(window* window, size_t kmer_size, phylo_kmer::score_type threshold) noexcept;
            DCLA_iterator(const DCLA_iterator&) = delete;
            DCLA_iterator(DCLA_iterator&&) = default;
            DCLA_iterator& operator=(const DCLA_iterator&) = delete;
            DCLA_iterator& operator=(DCLA_iterator&& rhs) noexcept;
            ~DCLA_iterator() noexcept override = default;

            DCLA_iterator& operator++() override;

        private:
            void preprocess();

            size_t _prefix_size;
        };

        /// Creates a phylo-k-mer enumerator for the window based on the algorithm selected
        std::unique_ptr<kmer_iterator_pimpl> make_iterator(const algorithm& algorithm,
                                                           window* window, size_t start_pos,
                                                           size_t kmer_size, phylo_kmer::score_type threshold,
                                                           impl::vector_type<unpositioned_phylo_kmer> prefixes);*/

    }

    /*
    /// Iterates over alive k-mers in a given window.
    /// This is just a composition-based wrapper for impl::kmer_iterator_pimpl
    class kmer_iterator
    {
    public:
        using reference = impl::kmer_iterator_pimpl::reference;
        using pointer = impl::kmer_iterator_pimpl::pointer;

        kmer_iterator(const xpas::algorithm& algorithm, window* window,
                      size_t kmer_size, phylo_kmer::score_type threshold);
        kmer_iterator(const kmer_iterator&) = delete;
        kmer_iterator(kmer_iterator&&) = default;
        kmer_iterator& operator=(const kmer_iterator&) = delete;
        kmer_iterator& operator=(kmer_iterator&& rhs) noexcept = default;
        ~kmer_iterator() noexcept = default;

        reference operator*() const noexcept;
        pointer operator->() const noexcept;

        bool operator==(const kmer_iterator& rhs) const noexcept;
        bool operator!=(const kmer_iterator& rhs) const noexcept;

        kmer_iterator& operator++();
    private:
        std::unique_ptr<impl::kmer_iterator_pimpl> _pimpl;
    };

    /// A range-based wrapper for kmer_iterator
    class enumerate_kmers
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using const_iterator = kmer_iterator;
        using reference = kmer_iterator::reference;

        enumerate_kmers(const algorithm& algorithm,
                        window* window, size_t kmer_size, phylo_kmer::score_type threshold,
                        impl::vector_type<unpositioned_phylo_kmer> prefixes);
        enumerate_kmers(const enumerate_kmers&) = delete;
        enumerate_kmers(enumerate_kmers&&) = delete;
        enumerate_kmers& operator=(const enumerate_kmers&) = delete;
        enumerate_kmers& operator=(enumerate_kmers&& other) = delete;
        ~enumerate_kmers() noexcept = default;

        [[nodiscard]]
        const_iterator begin() const;

        [[nodiscard]]
        const_iterator end() const noexcept;

    private:
        algorithm _algorithm;
        window* _window;
        size_t _kmer_size;
        phylo_kmer::score_type _threshold;
    };*/


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
}



#endif //RAPPAS2_PK_COMPUTE_H
