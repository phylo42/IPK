#ifndef XPAS_WINDOW_H
#define XPAS_WINDOW_H

#include "row.h"
#include "pk_compute.h"
#include <stack>
#include <memory>

namespace xpas::impl
{
    template<class T>
    using vector_type = std::vector<T>;
}

namespace xpas
{
    class node_entry;
    class window;

    namespace impl
    {

        /// An abstract class for phylo-k-mer iterators. Used as part of
        /// "pointer to implementation" idiom in kmer_iterator
        /// An iterator takes a window and implements an algorithm that
        /// enumerates all alive k-mers in this window given the threshold.
        /// Needed to be able to plug in and compare different phylo-k-mer
        /// computation algorithms
        class kmer_iterator_pimpl
        {
        public:
            /// Member types
            using iterator_category = std::forward_iterator_tag;
            using reference = const xpas::unpositioned_phylo_kmer&;
            using pointer = const xpas::unpositioned_phylo_kmer*;

            /// Constructors and a virtual destructor
            kmer_iterator_pimpl(window* window, size_t kmer_size, phylo_kmer::score_type threshold);
            virtual ~kmer_iterator_pimpl() noexcept = default;
            kmer_iterator_pimpl(const kmer_iterator_pimpl&) = delete;
            kmer_iterator_pimpl(kmer_iterator_pimpl&&) = default;

            /// Only move assignment is allowed
            kmer_iterator_pimpl& operator=(const kmer_iterator_pimpl&) = delete;
            /// not virtual
            kmer_iterator_pimpl& operator=(kmer_iterator_pimpl&& rhs) noexcept;

            bool operator==(const kmer_iterator_pimpl& rhs) const noexcept;
            bool operator!=(const kmer_iterator_pimpl& rhs) const noexcept;

            virtual kmer_iterator_pimpl& operator++() = 0;

            reference operator*() const noexcept;
            pointer operator->() const noexcept;

        protected:
            window* _window;

            size_t _kmer_size;
            phylo_kmer::score_type _threshold;

            xpas::unpositioned_phylo_kmer _current;
        };


        /// \brief Enumerates phylo-k-mers in a given window using divide-and-conquer
        /// (without applying the lookahead technique)
        class naive_DC_iterator : public kmer_iterator_pimpl
        {
        public:
            naive_DC_iterator(window* view, size_t kmer_size, phylo_kmer::score_type threshold,
                              phylo_kmer::pos_type start_pos, size_t prefix_size,
                              vector_type<unpositioned_phylo_kmer> prefixes) noexcept;
            naive_DC_iterator(const naive_DC_iterator&) = delete;
            naive_DC_iterator(naive_DC_iterator&&) = default;
            naive_DC_iterator& operator=(const naive_DC_iterator&) = delete;
            naive_DC_iterator& operator=(naive_DC_iterator&& rhs) noexcept;
            ~naive_DC_iterator() noexcept override = default;

            //bool operator==(const naive_dac_enumerator& rhs) const noexcept;
            //bool operator!=(const naive_dac_enumerator& rhs) const noexcept;

            naive_DC_iterator& operator++() override;
        private:
            unpositioned_phylo_kmer _next_kmer();

            void _select_suffix_bound();
            void _finish_iterator();

            size_t _prefix_size;
            phylo_kmer::pos_type _start_pos;

            vector_type<unpositioned_phylo_kmer> _prefixes;
            vector_type<unpositioned_phylo_kmer>::iterator _prefix_it;

            vector_type<unpositioned_phylo_kmer> _suffixes;
            vector_type<unpositioned_phylo_kmer>::iterator _suffix_it;
            vector_type<unpositioned_phylo_kmer>::iterator _last_suffix_it;
        };

        naive_DC_iterator make_dac_end_iterator();

        /// Creates a phylo-k-mer enumerator for the window based on the algorithm selected
        std::unique_ptr<kmer_iterator_pimpl> make_iterator(const algorithm& algorithm,
                                                           window* window, size_t start_pos,
                                                           size_t kmer_size, phylo_kmer::score_type threshold,
                                                           impl::vector_type<unpositioned_phylo_kmer> prefixes);

    }

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
    };

    /// \brief A lightweight view of node_entry. Implements a window of size k over a node_entry
    class window final
    {
    public:
        window(const node_entry* entry, phylo_kmer::score_type threshold,
               phylo_kmer::pos_type start, phylo_kmer::pos_type end) noexcept;
        window(const window& other) noexcept;
        window(window&&) = default;
        window& operator=(const window&) = delete;
        window& operator=(window&& other) noexcept;
        ~window() noexcept = default;

        [[nodiscard]]
        const node_entry* get_entry() const noexcept;

        [[nodiscard]]
        phylo_kmer::pos_type get_start_pos() const noexcept;
        void set_start_pos(phylo_kmer::pos_type pos);

        [[nodiscard]]
        phylo_kmer::pos_type get_end_pos() const noexcept;
        void set_end_pos(phylo_kmer::pos_type pos);

        [[nodiscard]]
        phylo_kmer::score_type get_threshold() const noexcept;

        void set_prefixes(impl::vector_type<unpositioned_phylo_kmer> prefixes);

        [[nodiscard]]
        size_t get_prefix_size() const noexcept;

        void set_prefix_size(size_t size);

    private:
        const node_entry* _entry;

        phylo_kmer::score_type _threshold;
        phylo_kmer::pos_type _start;
        phylo_kmer::pos_type _end;

        /// The vector of precomputed prefixes for this window.
        /// Can be obtained from iteration of the previous window
        /// (according to the order of windows, given by chain_windows)
        impl::vector_type<xpas::unpositioned_phylo_kmer> _prefixes;

        /// The length of the precomputed prefixes stored in _prefixes.
        /// It is zero if _prefixes if empty.
        size_t _prefix_size;
    };

    bool operator==(const xpas::window& a, const xpas::window& b) noexcept;
    bool operator!=(const xpas::window& a, const xpas::window& b) noexcept;

}

#endif //XPAS_WINDOW_H
