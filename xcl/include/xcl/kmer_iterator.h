#ifndef XPAS_KMER_ITERATOR_H
#define XPAS_KMER_ITERATOR_H

#include <string_view>
#include <vector>
#include <tuple>
#include "phylo_kmer.h"

namespace xpas
{
    /// \brief An iterator class for k-mers of an input sequence.
    /// \details Iterates over all the k-mers of a given sequence and encodes them
    /// in a rolling fashion. Returns a lightweight std::string_view to a current k-mer
    /// of a given sequence and its code.
    template<typename AmbiguityPolicy>
    class kmer_iterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;

        /// Return value of a kmer_iterator, which is a pair of a string view to a k-mer
        /// and a k-mer code storage
        using value_type = std::pair<std::string_view, typename AmbiguityPolicy::value_type>;

        using pos_t = size_t;

        kmer_iterator(std::string_view sequence_view, size_t kmer_size) noexcept;
        kmer_iterator(const kmer_iterator&) = delete;
        kmer_iterator(kmer_iterator&&) = delete;
        kmer_iterator& operator=(const kmer_iterator&) = delete;
        kmer_iterator& operator=(kmer_iterator&&) = delete;
        ~kmer_iterator() noexcept = default;

        kmer_iterator& operator++()
        {
            bool last_kmer_is_valid = true;
            bool stop = false;
            while (!stop)
            {
                /// if the last k-mer is valid, we can reuse its code in a rolling fashion
                /// otherwise we need to calculate the code from scratch
                const auto encode_result = last_kmer_is_valid ? encode_kmer_from_previous() : encode_kmer();
                std::tie(_kmer_value, _kmer_position, last_kmer_is_valid) = encode_result;
                stop = last_kmer_is_valid || (_kmer_position == std::string_view::npos);
            }
            return *this;
        }

        bool operator==(const kmer_iterator& rhs) const noexcept
        {
            return (_kmer_position == rhs._kmer_position) && (_sequence_view == rhs._sequence_view);
        }

        bool operator!=(const kmer_iterator& rhs) const noexcept
        {
            return !(*this == rhs);
        }

        value_type operator*() const noexcept
        {
            return { _sequence_view.substr(_kmer_position, _kmer_size), _kmer_value };
        }

    private:
        void init()
        {
            /// _kmer_size == 0 means npos
            if (_kmer_size == 0)
            {
                _kmer_position = std::string_view::npos;
            }
            else
            {
                /// calculate the code of the first valid k-mer
                bool stop = false;
                while (!stop)
                {
                    bool is_valid = false;
                    std::tie(_kmer_value, _kmer_position, is_valid) = encode_kmer();
                    stop = is_valid || (_kmer_position == std::string_view::npos);
                }
            }
        }

        std::tuple<typename AmbiguityPolicy::value_type, pos_t, bool> encode_kmer();
        std::tuple<typename AmbiguityPolicy::value_type, pos_t, bool> encode_kmer_from_previous() const;

        /// String view of the original sequence
        std::string_view _sequence_view;
        size_t _kmer_size;

        /// A current stored value of iterator, which type depends on the policy
        /// (Ex. one or more k-mer codes)
        typename AmbiguityPolicy::value_type _kmer_value;

        /// A current position of the iterator in the input sequence
        pos_t _kmer_position;
    };

    /// \brief A little proxy class that just creates a k-mer iterator over a sequence.
    /// Example: for (const auto& [kmer, code] : to_kmers<no_ambiguity_policy>(sequence_view, kmer_size)) { ... }
    template<typename AmbiguityPolicy>
    class to_kmers
    {
    public:
        using const_iterator = kmer_iterator<AmbiguityPolicy>;

        to_kmers(std::string_view sequence_view, size_t kmer_size) noexcept
            : _sequence_view{ sequence_view }, _kmer_size{ kmer_size }
        {
            /// If the input string is too small for a given k, we force begin() == end()
            /// to represent an empty set of k-mers.
            if (_sequence_view.size() < _kmer_size)
            {
                _kmer_size = 0;
            }
        }

        const_iterator begin() const
        {
            return { _sequence_view, _kmer_size };
        }

        const_iterator end() const
        {
            return { _sequence_view, 0 };
        }

    private:
        std::string_view _sequence_view;
        size_t _kmer_size;
    };
}

#endif //RAPPAS_CORE_KMER_ITERATOR_H
