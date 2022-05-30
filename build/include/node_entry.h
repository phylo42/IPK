#ifndef RAPPAS_CPP_NODE_ENTRY_H
#define RAPPAS_CPP_NODE_ENTRY_H

#include <vector>
#include "window.h"

namespace xpas
{
    class view_iterator;

    /// \brief A submatrix of posterior probabilities matrix (fixed branch, all the positions of input alignment)
    class node_entry final
    {
    public:
        using vector_type = std::vector<xpas::row_type>;

        explicit node_entry() noexcept = default;
        node_entry(std::string _id, vector_type&& rows);
        node_entry(const node_entry&) = delete;
        node_entry(node_entry&&) = default;
        node_entry& operator=(const node_entry&) = delete;
        node_entry& operator=(node_entry&&) = default;
        ~node_entry() noexcept = default;

        void push_back(xpas::row_type row);

        [[nodiscard]]
        size_t get_alignment_size() const;

        [[nodiscard]]
        std::string get_label() const;

        [[nodiscard]]
        const xpas::proba_pair& at(size_t position, size_t variant) const;

    private:
        std::string _branch_label;
        vector_type _rows;
    };

    namespace impl
    {
        /// Chains windows of the posterior probability matrix in multiple chains.
        /// For two consecutive windows of a chain, the suffixes of the predecessor are
        /// used as prefixes by the successor (see dac_kmer_iterator).
        /// The distance between two windows is floor(k/2) or ceil(k/2) depending on
        /// the parity of the number of the window in the chain.
        ///
        /// Thus, it chains windows in the following manner:
        /// 0) [0 .. k-1] -> [floor(k/2) .. k+floor(k/2)-1] -> [k .. 2k-1] -> [k+floor(k/2) .. 2k+floor(k/2)-1] -> ...
        /// ...
        /// i) [i .. i+k-1] -> [i+floor(k/2) .. i+k+floor(k/2)-1] -> [i+k .. i+2k-1] -> ...
        /// ...
        /// The last value of i is floor(k/2)-1, which means there are floor(k/2) chains.
        ///
        class chain_window_iterator
        {
        public:
            using iterator_category = std::forward_iterator_tag;
            using reference = window&;

            chain_window_iterator(window view, size_t kmer_size, phylo_kmer::score_type threshold) noexcept;
            chain_window_iterator(const chain_window_iterator&) = delete;
            chain_window_iterator(chain_window_iterator&&) = delete;
            chain_window_iterator& operator=(const chain_window_iterator&) = delete;
            chain_window_iterator& operator=(chain_window_iterator&&) = delete;
            ~chain_window_iterator() = default;

            chain_window_iterator& operator++();

            bool operator==(const chain_window_iterator& rhs) const noexcept;
            bool operator!=(const chain_window_iterator& rhs) const noexcept;

            reference operator*() noexcept;
        private:
            window _view;

            size_t _kmer_size;
            phylo_kmer::score_type _threshold;

            /// The position of the first windows in the current chain,
            /// which is "i" in the notation above
            size_t _first_view_pos;
        };
    }

    /// Constructing iterator of node_entry
    class chain_windows
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using const_iterator = impl::chain_window_iterator;

        using reference = window&;

        chain_windows(const node_entry& entry, size_t kmer_size, phylo_kmer::score_type threshold);
        chain_windows(const chain_windows&) = delete;
        chain_windows(chain_windows&&) = delete;
        chain_windows& operator=(const chain_windows&) = delete;
        chain_windows& operator=(chain_windows&&) = delete;
        ~chain_windows() noexcept = default;

        [[nodiscard]]
        const_iterator begin() const;

        [[nodiscard]]
        const_iterator end() const noexcept;
    private:
        const node_entry& _entry;

        size_t _kmer_size;
        xpas::phylo_kmer::pos_type _start_pos;
        xpas::phylo_kmer::score_type _threshold;
    };

    bool operator==(const node_entry& lhs, const node_entry& rhs);
}


#endif //RAPPAS_CPP_NODE_ENTRY_H
