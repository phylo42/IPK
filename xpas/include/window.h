#ifndef XPAS_WINDOW_H
#define XPAS_WINDOW_H

#include "row.h"
#include <stack>
#include <memory>

namespace xpas::impl
{
    template<class T>
    using vector_type = std::vector<T>;

    using score_t = phylo_kmer::score_type;
}

namespace xpas
{
    class matrix
    {
    public:
        using column = impl::vector_type<impl::score_t>;

        matrix() noexcept = default;
        matrix(std::vector<column> data, std::string label);
        matrix(const matrix&) = delete;
        matrix(matrix&&) noexcept = default;
        ~matrix() noexcept = default;

        matrix& operator=(const matrix&) = delete;
        matrix& operator=(matrix&&) noexcept = default;

        /// Range product query for the max values of every column
        void preprocess();

        [[nodiscard]]
        impl::score_t get(size_t i, size_t j) const;

        [[nodiscard]]
        size_t width() const;

        [[nodiscard]]
        bool empty() const;

        void set_label(const std::string& label);

        [[nodiscard]]
        const std::vector<column>& get_data() const;

        std::vector<column>& get_data();

        [[nodiscard]]
        std::string get_label() const;

        [[nodiscard]]
        const column& get_column(size_t j) const;

        [[nodiscard]]
        impl::score_t range_max_sum(size_t start_pos, size_t len) const;

    private:
        std::vector<column> _data;

        std::string _label;

        std::vector<impl::score_t> _best_scores;
    };

    namespace impl
    {
        class window_iterator;
        class chained_window_iterator;
    }

    class window
    {
        friend class impl::window_iterator;
        friend class impl::chained_window_iterator;
    public:
        window(const matrix* m, size_t start_pos, size_t size);
        window(const window&) = delete;
        window(window&&) noexcept = default;
        //window& operator=(const window& other) = default;
        window& operator=(const window& other) = delete;
        window& operator=(window&&) noexcept = default;

        bool operator==(const window& other) const;
        bool operator!=(const window& other) const;

        [[nodiscard]]
        impl::score_t get(size_t i, size_t j) const;

        [[nodiscard]]
        size_t size() const;

        [[nodiscard]]
        bool empty() const;

        [[nodiscard]]
        size_t get_position() const;

        [[nodiscard]]
        phylo_kmer::score_type range_max_product(size_t pos, size_t len) const;

        [[nodiscard]]
        std::pair<size_t, impl::score_t> max_at(size_t column) const;

        [[nodiscard]]
        matrix::column get_column(size_t j) const;


    private:
        const matrix* _matrix;
        size_t _start_pos;
        size_t _size;
    };

    namespace impl
    {
        class window_iterator
        {
        public:
            using iterator_category = std::forward_iterator_tag;
            using reference = window&;

            window_iterator(const matrix* matrix, size_t kmer_size) noexcept;
            window_iterator(const window_iterator&) = delete;
            window_iterator(window_iterator&&) = delete;
            window_iterator& operator=(const window_iterator&) = delete;
            window_iterator& operator=(window_iterator&&) = delete;
            ~window_iterator() = default;

            window_iterator& operator++();

            bool operator==(const window_iterator& rhs) const noexcept;
            bool operator!=(const window_iterator& rhs) const noexcept;

            reference operator*() noexcept;
        private:
            const matrix* _matrix;

            window _window;

            size_t _kmer_size;

            size_t _current_pos;
        };

        class chained_window_iterator
        {
        public:
            using iterator_category = std::forward_iterator_tag;
            using reference = window&;

            chained_window_iterator(const matrix* matrix, size_t kmer_size);
            chained_window_iterator(const chained_window_iterator&) = delete;
            chained_window_iterator(window_iterator&&) = delete;
            chained_window_iterator& operator=(const chained_window_iterator&) = delete;
            chained_window_iterator& operator=(chained_window_iterator&&) = delete;
            ~chained_window_iterator() = default;

            chained_window_iterator& operator++();

            bool operator==(const chained_window_iterator& rhs) const noexcept;
            bool operator!=(const chained_window_iterator& rhs) const noexcept;

            std::tuple<reference, reference, reference> operator*() noexcept;
        private:
            window _get_next_window();

            const matrix* _matrix;

            window _window;
            window _previous_window;
            window _next_window;

            size_t _kmer_size;

            // the first position j of the current chain of windows
            size_t _chain_start;

            // the first position j of the last chain possible
            size_t _last_chain_pos;
        };
    }

    class to_windows
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using const_iterator = impl::window_iterator;

        using reference = window&;

        to_windows(const matrix* matrix, size_t kmer_size);
        to_windows(const to_windows&) = delete;
        to_windows(to_windows&&) = delete;
        to_windows& operator=(const to_windows&) = delete;
        to_windows& operator=(to_windows&&) = delete;
        ~to_windows() noexcept = default;

        [[nodiscard]]
        const_iterator begin() const;

        [[nodiscard]]
        const_iterator end() const noexcept;

    private:
        const matrix* _matrix;
        size_t _kmer_size;
    };

    class chain_windows
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using const_iterator = impl::chained_window_iterator;

        using reference = window&;

        chain_windows(const matrix* matrix, size_t kmer_size);
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
        const matrix* _matrix;
        size_t _kmer_size;
    };
}
#endif //XPAS_WINDOW_H
