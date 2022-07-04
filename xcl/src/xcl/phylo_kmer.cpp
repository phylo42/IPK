#include "xcl/phylo_kmer.h"
#include "xcl/seq.h"
#include <cmath>
#include <vector>

using namespace xcl;

/// Compares two floats for almost equality.
/// From: https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type almost_equal(T x, T y, int ulp = 1) noexcept
{
    // The machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x - y) <= std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
           // unless the result is subnormal
           || std::abs(x - y) < std::numeric_limits<T>::min();
}

template <>
positioned_pkdb_value xcl::make_pkdb_value<positioned_phylo_kmer>(phylo_kmer::branch_type branch,
                                                                   phylo_kmer::score_type score,
                                                                   phylo_kmer::pos_type position)
{
    return { branch, score, position };
}

template <>
unpositioned_pkdb_value xcl::make_pkdb_value<unpositioned_phylo_kmer>(phylo_kmer::branch_type branch,
                                                                       phylo_kmer::score_type score,
                                                                       phylo_kmer::pos_type position)
{
    (void)position;
    return { branch, score };
}

template <>
positioned_phylo_kmer xcl::make_phylo_kmer(phylo_kmer::key_type key, phylo_kmer::score_type score,
                                            phylo_kmer::pos_type position)
{
    return { key, score, position };
}

template <>
unpositioned_phylo_kmer xcl::make_phylo_kmer(phylo_kmer::key_type key, phylo_kmer::score_type score,
                                              phylo_kmer::pos_type position)
{
    (void)position;
    return { key, score };
}

bool xcl::operator==(const positioned_phylo_kmer& lhs, const positioned_phylo_kmer& rhs) noexcept
{
    if (lhs.is_nan() || rhs.is_nan())
    {
        return false;
    }
    else
    {
        return (lhs.key == rhs.key) && (almost_equal<phylo_kmer::score_type>(lhs.score, rhs.score)) &&
            lhs.position == rhs.position;
    }
}

bool xcl::operator==(const unpositioned_phylo_kmer& lhs, const unpositioned_phylo_kmer& rhs) noexcept
{
    if (lhs.is_nan() || rhs.is_nan())
    {
        return false;
    }
    else
    {
        return (lhs.key == rhs.key) && (almost_equal<phylo_kmer::score_type>(lhs.score, rhs.score));
    }
}

phylo_kmer::score_type xcl::score_threshold(phylo_kmer::score_type omega, size_t kmer_size)
{
    return std::pow(omega / seq_traits::alphabet_size, phylo_kmer::score_type(kmer_size));
}

template<>
optional<no_ambiguity_policy::value_type> xcl::encode_kmer<no_ambiguity_policy>(std::string_view kmer)
{
    no_ambiguity_policy::value_type key = 0;
    for (const auto base : kmer)
    {
        if (const auto& base_code = encode<no_ambiguity_policy>(base); base_code)
        {
            key <<= bit_length<seq_type>();
            key |= *base_code;
        }
        else
        {
            return nullopt;
        }
    }
    return key;
}

template<>
optional<one_ambiguity_policy::value_type> xcl::encode_kmer<one_ambiguity_policy>(std::string_view kmer)
{
    auto keys = one_ambiguity_policy::value_type();
    keys.push_back(0);

    size_t num_ambiguities = 0;

    for (const auto base : kmer)
    {
        if (const auto& base_codes = encode<one_ambiguity_policy>(base); base_codes)
        {
            /// if the character is ambiguous
            if (base_codes->size() > 1)
            {
                /// allow only one ambiguity per k-mer
                if (num_ambiguities > 0)
                {
                    return nullopt;
                }
                else
                {
                    auto old_key = keys.back();
                    keys.pop_back();
                    for (const auto base_code : *base_codes)
                    {
                        auto new_key = old_key;
                        new_key <<= bit_length<seq_type>();
                        new_key |= base_code;
                        keys.push_back(new_key);
                    }
                }
                ++num_ambiguities;
            }
            /// if not ambiguous
            else if (base_codes->size() == 1)
            {
                for (auto& key : keys)
                {
                    key <<= bit_length<seq_type>();
                    key |= (*base_codes)[0];
                }
            }
        }
        else
        {
            return nullopt;
        }
    }
    return keys;
}


std::string xcl::decode_kmer(phylo_kmer::key_type key, size_t kmer_size)
{
    std::vector<uint8_t> result;
    result.reserve(kmer_size);

    while (key > 0)
    {
        result.push_back(decode(key & ~rightmost_symbol_mask<seq_type>()));
        key >>= bit_length<seq_type>();
    }

    for (size_t i = result.size(); i < kmer_size; ++i)
    {
        result.push_back(decode(0));
    }

    return { result.rbegin(), result.rend() };
}

