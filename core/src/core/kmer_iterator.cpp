#include <tuple>
#include "core/kmer_iterator.h"
#include "core/phylo_kmer.h"

using namespace core;


template<>
std::tuple<no_ambiguity_policy::value_type,
    kmer_iterator<no_ambiguity_policy>::pos_t,
    bool> kmer_iterator<no_ambiguity_policy>::encode_kmer()
{
    phylo_kmer::key_type key = phylo_kmer::nan_key;
    bool success = false;

    /// string is too small to start from _kmer_position
    if (_kmer_position > _sequence_view.size() - _kmer_size)
    {
        _kmer_position = std::string_view::npos;
    }
    else
    {
        const auto kmer_view = _sequence_view.substr(_kmer_position, _kmer_size);

        /// correct k-mer with gaps
        if (auto code = core::encode_kmer(kmer_view); code)
        {
            key = *code;
            success = true;
        }
            /// skip a k-mer if we found a gap or ambiguity
        else
        {
            /// WARNING:
            /// this is not effective, because we can skip all the bases up to the one
            /// that is not valid, not just only one. But this is not clear how to return
            /// this position from core::encode_kmer, still having the API of the core clear and simple.
            /// TODO: find a solution
            _kmer_position += 1;
        }
    }
    return { key, _kmer_position, success };
}

template <>
kmer_iterator<no_ambiguity_policy>::kmer_iterator(std::string_view sequence_view, size_t kmer_size) noexcept
    : _sequence_view{ sequence_view }, _kmer_size{ kmer_size }, _kmer_value{ phylo_kmer::nan_key }, _kmer_position{ 0 }
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

template<>
kmer_iterator<no_ambiguity_policy>::value_type kmer_iterator<no_ambiguity_policy>::operator*() const noexcept
{
    return { _sequence_view.substr(_kmer_position, _kmer_size), _kmer_value };
}

template<>
std::tuple<no_ambiguity_policy::value_type,
    kmer_iterator<no_ambiguity_policy>::pos_t,
    bool> kmer_iterator<no_ambiguity_policy>::encode_kmer_from_previous() const
{
    auto key = _kmer_value;
    auto position = _kmer_position;
    bool success = false;

    if (position <= _sequence_view.size() - _kmer_size)
    {
        /// get the last base of the new k-mer
        ++position;
        const auto new_base = _sequence_view[position + _kmer_size - 1];

        /// check for non-valid symbols (gaps, '*', '.' etc.)
        if (auto new_code = core::encode(new_base); new_code)
        {
            /// shift bits and remove the first base of the current k-mer
            key <<= core::bit_length<core::seq_type>();
            key &= ((1u << (core::bit_length<core::seq_type>() * _kmer_size)) - 1);

            /// add the new base
            key |= *new_code;
            success = true;
        }
            /// if we found a gap or ambiguity, skip it and go to the next possible k-mer
        else
        {
            position = position + _kmer_size;
        }
    }
        /// if we reached the end of sequence
    else
    {
        position = std::string_view::npos;
    }
    return { key, position, success };
}


to_kmers::to_kmers(std::string_view sequence_view, size_t kmer_size) noexcept
    : _sequence_view{ sequence_view }, _kmer_size{ kmer_size }
{
    /// If the input string is too small for a given k, we force begin() == end()
    /// to represent an empty set of k-mers.
    if (_sequence_view.size() < _kmer_size)
    {
        _kmer_size = 0;
    }
}

to_kmers::const_iterator to_kmers::begin() const
{
    return { _sequence_view, _kmer_size };
}

to_kmers::const_iterator to_kmers::end() const
{
    return { _sequence_view, 0 };
}