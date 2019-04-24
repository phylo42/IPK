#include "core/kmer_iterator.h"
#include "core/phylo_kmer.h"

using namespace core;

kmer_iterator::kmer_iterator(std::string_view sequence_view, size_t kmer_size) noexcept
    : _sequence_view{ sequence_view }
    , _kmer_position{ 0 }
    , _kmer_size{ kmer_size }
{
    /// _kmer_size == 0 means npos
    if (_kmer_size == 0)
    {
        _kmer_position = std::string_view::npos;
    }
    else
    {
        _current_key = core::encode_kmer(_sequence_view.substr(_kmer_position, _kmer_size));
    }
}

kmer_iterator& kmer_iterator::operator++()
{
    if (_kmer_position < _sequence_view.size() - _kmer_size)
    {
        /// shift bits and remove the first base of the current k-mer
        _current_key <<= core::bit_length<core::seq_type>();
        _current_key &= ((1u << (core::bit_length<core::seq_type>() * _kmer_size)) - 1);

        // add the last base of the new k-mer
        ++_kmer_position;
        const auto new_base = _sequence_view[_kmer_position + _kmer_size - 1];
        _current_key |= core::encode(new_base);
    }
    else
    {
        _kmer_position = std::string_view::npos;
    }
    return *this;
}

bool kmer_iterator::operator==(const kmer_iterator& rhs) const noexcept
{
    return (_kmer_position == rhs._kmer_position) && (_sequence_view == rhs._sequence_view);
}

bool kmer_iterator::operator!=(const kmer_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

kmer_iterator::value_type kmer_iterator::operator*() const noexcept
{
    return { _sequence_view.substr(_kmer_position, _kmer_size), _current_key };
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