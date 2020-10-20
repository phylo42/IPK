#include <tuple>
#include "xpas/kmer_iterator.h"
#include "xpas/phylo_kmer.h"

using namespace xpas;

/// No ambiguity policy
template<>
std::tuple<no_ambiguity_policy::value_type,
    kmer_iterator<no_ambiguity_policy>::pos_t,
    bool> kmer_iterator<no_ambiguity_policy>::encode_kmer()
{
    phylo_kmer::key_type key = phylo_kmer::na_key;
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
        if (auto code = xpas::encode_kmer<no_ambiguity_policy>(kmer_view); code)
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
    : _sequence_view{ sequence_view }, _kmer_size{ kmer_size }, _kmer_value{ phylo_kmer::na_key }, _kmer_position{ 0 }
{
   init();
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
        if (auto new_code = xpas::encode<no_ambiguity_policy>(new_base); new_code)
        {
            /// shift bits and remove the first base of the current k-mer
            key <<= xpas::bit_length<xpas::seq_type>();
            key &= ((1u << (xpas::bit_length<xpas::seq_type>() * _kmer_size)) - 1);

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

/// One ambiguity policy
template<>
std::tuple<one_ambiguity_policy::value_type,
    kmer_iterator<one_ambiguity_policy>::pos_t,
    bool> kmer_iterator<one_ambiguity_policy>::encode_kmer()
{
    one_ambiguity_policy::value_type keys;
    bool success = false;

    /// string is too small to start from _kmer_position
    if (_kmer_position > _sequence_view.size() - _kmer_size)
    {
        _kmer_position = std::string_view::npos;
    }
    else
    {
        const auto kmer_view = _sequence_view.substr(_kmer_position, _kmer_size);

        /// resolve ambiguities
        if (auto codes = xpas::encode_kmer<one_ambiguity_policy>(kmer_view); codes)
        {
            keys = *codes;
            success = true;
        }
        /// skip a k-mer if we found a gap or an unrelated character
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
    return { keys, _kmer_position, success };
}

template<>
std::tuple<one_ambiguity_policy::value_type,
    kmer_iterator<one_ambiguity_policy>::pos_t,
    bool> kmer_iterator<one_ambiguity_policy>::encode_kmer_from_previous() const
{
    auto keys = _kmer_value;
    auto position = _kmer_position;
    bool success = false;

    /// If not the end of sequence
    if (position <= _sequence_view.size() - _kmer_size)
    {
        /// get the last base of the new k-mer
        ++position;
        const auto new_base = _sequence_view[position + _kmer_size - 1];

        /// encode a new base resolving ambiguities
        auto new_codes = xpas::encode<one_ambiguity_policy>(new_base);

        /// Calculate if the (k-1) suffix of the previous k-mer is still ambiguous
        bool is_ambiguous = false;
        if (keys.size() > 1)
        {
            auto suffix1 = keys[0];
            suffix1 <<= xpas::bit_length<xpas::seq_type>();
            suffix1 &= ((1u << (xpas::bit_length<xpas::seq_type>() * _kmer_size)) - 1);

            auto suffix2 = keys[1];
            suffix2 <<= xpas::bit_length<xpas::seq_type>();
            suffix2 &= ((1u << (xpas::bit_length<xpas::seq_type>() * _kmer_size)) - 1);

            is_ambiguous = (suffix1 != suffix2);
        }

        /// check if is a second ambiguity in the k-mer
        bool more_than_1_ambiguity = new_codes && is_ambiguous && (new_codes->size() > 1);
        if (new_codes && !more_than_1_ambiguity)
        {
            if (new_codes->size() == 1)
            {
                for (auto& key : keys)
                {
                    for (const auto& new_code : *new_codes)
                    {
                        /// shift bits and remove the first base of the current k-mer
                        key <<= xpas::bit_length<xpas::seq_type>();
                        key &= ((1u << (xpas::bit_length<xpas::seq_type>() * _kmer_size)) - 1);

                        /// add the new base
                        key |= new_code;
                    }
                }

                /// If there was an ambiguity at the beginning that we cut off
                /// without adding an ambiguous base at the end, the new k-mer
                /// is not ambiguous, so we need to shrink down the vector of codes
                /// Example: the previous k-mer was .AAA, we had 4 codes
                /// the new one is AAAT, we have only one code
                if (keys.size() > 1 && keys[0] == keys[1])
                {
                    keys.erase(keys.begin() + 1, keys.end());
                }
            }
            else
            {
                auto key = keys.back();
                /// all the old k-mers have the same suffix here. We need only one of them.
                keys.clear();

                for (const auto new_code : *new_codes)
                {
                    /// shift bits and remove the first base of the current k-mer
                    auto new_key = key;
                    new_key <<= xpas::bit_length<xpas::seq_type>();
                    new_key &= ((1u << (xpas::bit_length<xpas::seq_type>() * _kmer_size)) - 1);

                    /// add the new base
                    new_key |= new_code;
                    keys.push_back(new_key);
                }
            }
            success = true;
        }
        /// If we found a second ambiguity, skip it a position. Here we can not skip
        /// kmer_size positions, because the next k-mer may have just only one ambiguity
        else
        {
            position = position + 1;
        }
    }
    /// if we reached the end of sequence
    else
    {
        position = std::string_view::npos;
    }
    return { keys, position, success };
}

template <>
kmer_iterator<one_ambiguity_policy>::kmer_iterator(std::string_view sequence_view, size_t kmer_size) noexcept
    : _sequence_view{ sequence_view }, _kmer_size{ kmer_size }, _kmer_value{ }, _kmer_position{ 0 }
{
    init();
}