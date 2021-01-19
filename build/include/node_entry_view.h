#ifndef XPAS_NODE_ENTRY_VIEW_H
#define XPAS_NODE_ENTRY_VIEW_H

#include "row.h"

namespace xpas::impl
{
    template<class T>
    using vector_type = std::vector<T>;
}

namespace xpas
{
    class node_entry;
    class node_entry_view;

    namespace impl
    {
        /// \brief Temporary data storage for mmers, m <= k. Used only in the branch-and-bound algorithm
        /// to generate phylo-kmers
        struct phylo_mmer
        {
            xpas::unpositioned_phylo_kmer mmer;
            /// a position of the last letter
            xpas::phylo_kmer::pos_type last_position;
            size_t last_index;
            size_t next_index;
        };

        /// \brief Divide-and-conquer phylo-kmer iterator.
        class dac_kmer_iterator
        {
        public:
            /// Member types
            using iterator_category = std::forward_iterator_tag;
            using reference = const xpas::unpositioned_phylo_kmer&;
            using pointer = const xpas::unpositioned_phylo_kmer*;

            dac_kmer_iterator(node_entry_view* view, size_t kmer_size, xpas::phylo_kmer::score_type threshold,
                              xpas::phylo_kmer::pos_type start_pos, size_t prefix_size,
                              vector_type<xpas::unpositioned_phylo_kmer> prefixes) noexcept;
            dac_kmer_iterator(const dac_kmer_iterator&) = delete;
            dac_kmer_iterator(dac_kmer_iterator&&) = default;
            dac_kmer_iterator& operator=(const dac_kmer_iterator&) = delete;
            dac_kmer_iterator& operator=(dac_kmer_iterator&& rhs) noexcept;
            ~dac_kmer_iterator() noexcept = default;

            bool operator==(const dac_kmer_iterator& rhs) const noexcept;
            bool operator!=(const dac_kmer_iterator& rhs) const noexcept;
            dac_kmer_iterator& operator++();

            reference operator*() const noexcept;
            pointer operator->() const noexcept;
        private:
            xpas::unpositioned_phylo_kmer _next_phylokmer();
            void _select_suffix_bound();

            node_entry_view* _entry_view;

            size_t _kmer_size;
            size_t _prefix_size;
            xpas::phylo_kmer::pos_type _start_pos;
            xpas::phylo_kmer::score_type _threshold;
            xpas::unpositioned_phylo_kmer _current;

            vector_type<xpas::unpositioned_phylo_kmer> _prefixes;
            vector_type<xpas::unpositioned_phylo_kmer>::iterator _prefix_it;

            vector_type<xpas::unpositioned_phylo_kmer> _suffixes;
            vector_type<xpas::unpositioned_phylo_kmer>::iterator _suffix_it;
            vector_type<xpas::unpositioned_phylo_kmer>::iterator _last_suffix_it;
        };

        dac_kmer_iterator make_dac_end_iterator();
    }

    /// \brief A lightweight view of node_entry. Implements a "window" of size K over a node_entry.
    class node_entry_view final
    {
    public:
        using iterator = xpas::impl::dac_kmer_iterator;
        using reference = iterator::reference;

        node_entry_view(const node_entry* entry, xpas::phylo_kmer::score_type threshold,
                        xpas::phylo_kmer::pos_type start, xpas::phylo_kmer::pos_type end) noexcept;
        node_entry_view(const node_entry_view& other) noexcept;
        node_entry_view(node_entry_view&&) = default;
        node_entry_view& operator=(const node_entry_view&) = delete;
        node_entry_view& operator=(node_entry_view&& other) noexcept;
        ~node_entry_view() noexcept = default;

        [[nodiscard]]
        iterator begin();

        [[nodiscard]]
        iterator end() const noexcept;

        [[nodiscard]]
        const node_entry* get_entry() const noexcept;

        [[nodiscard]]
        phylo_kmer::pos_type get_start_pos() const noexcept;
        void set_start_pos(xpas::phylo_kmer::pos_type pos);

        [[nodiscard]]
        phylo_kmer::pos_type get_end_pos() const noexcept;
        void set_end_pos(phylo_kmer::pos_type pos);

        [[nodiscard]]
        phylo_kmer::score_type get_threshold() const noexcept;

        void set_prefixes(impl::vector_type<xpas::unpositioned_phylo_kmer> prefixes);

        [[nodiscard]]
        size_t get_prefix_size() const noexcept;

        void set_prefix_size(size_t size);

    private:
        const node_entry* _entry;

        xpas::phylo_kmer::score_type _threshold;
        xpas::phylo_kmer::pos_type _start;
        xpas::phylo_kmer::pos_type _end;

        /// The vector of precomputed prefixes for this window.
        /// Can be obtained from iteration of the previous window
        /// (according to the order of windows, given by chain_windows)
        impl::vector_type<xpas::unpositioned_phylo_kmer> _prefixes;

        /// The length of the precomputed prefixes stored in _prefixes.
        /// It is zero if _prefixes if empty.
        size_t _prefix_size;
    };

    bool operator==(const xpas::node_entry_view& a, const xpas::node_entry_view& b) noexcept;
    bool operator!=(const xpas::node_entry_view& a, const xpas::node_entry_view& b) noexcept;

}

#endif //XPAS_NODE_ENTRY_VIEW_H
