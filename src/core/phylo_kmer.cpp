#include "core/phylo_kmer.h"
#include <cmath>

using namespace core;

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

const phylo_kmer::key_type nan_key = 0;
const phylo_kmer::score_type nan_score = std::numeric_limits<phylo_kmer::score_type>::quiet_NaN();

bool phylo_kmer::is_nan() const
{
    return (key == nan_key) && (score == nan_score);
}

bool core::operator==(const phylo_kmer& lhs, const phylo_kmer& rhs) noexcept
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

phylo_kmer core::make_napk()
{
    return phylo_kmer { nan_key, nan_score };
}

phylo_kmer::score_type core::score_threshold(size_t kmer_size)
{
    return std::log10(powf(1.0f / seq_traits::alphabet_size, phylo_kmer::score_type(kmer_size)));
}

