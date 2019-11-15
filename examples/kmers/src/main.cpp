#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>
#include <core/seq.h>
#include <core/phylo_kmer.h>
#include <core/kmer_iterator.h>

/// Iterate over all the combinations with repetition
/// http://shoaib-ahmed.com/2018/for-each-combination-with-repetetion-c++/
template<typename V, typename Callable>
void for_each_combination(V &v, size_t gp_sz, Callable f) {
    V gp(gp_sz);
    auto total_n = std::pow(v.size(), gp.size());
    for (auto i = 0; i < total_n; ++i) {
        auto n = i;
        for (auto j = 0ul; j < gp.size(); ++j) {
            gp[gp.size() - j - 1] = v[n % v.size()];
            n /= v.size();
        }
        f(gp);
    }
}

void encode_unambiguous_string()
{
    //auto alphabet = std::vector<char>{ 'A', 'C', 'G', 'T', '-' };
    auto alphabet = std::vector<char>{ 'A', 'C', 'G', 'T'};
    const size_t kmer_size = 4;

    for_each_combination(alphabet, kmer_size,
                         [&](std::vector<char>& bases) {
                             const auto kmer = std::string{ bases.begin(), bases.end() };
                             if (const auto key = core::encode_kmer<core::no_ambiguity_policy>(kmer); key)
                             {
                                 std::cout << kmer << ": " << *key << std::endl;
                                 assert(kmer == core::decode_kmer(*key, kmer.size()));
                             }
                             else
                             {
                                 std::cout << kmer << " skipped" << std::endl;
                             }
                         });
}

void encode_unambiguous_string_view()
{
    std::cout << "\nIteration: " << std::endl;

    const auto long_read = std::string{ "--TTTAT-AAATGNNNN-CAAAN.NNTTTT---" };
    const size_t kmer_size = 4;

    /// core::to_kmers wrapper
    const auto to_kmers = [](std::string_view read, size_t k) { return core::to_kmers<core::no_ambiguity_policy>(read, k); };

    /// Iterate over k-mers of a string_view (implicitly created from long_read).
    /// This is more effective than calling core::encode_kmer for every k-mer of a string,
    /// because core::to_kmers calculates codes in a rolling fashion, if possible.
    for (const auto& [kmer, code] : to_kmers(long_read, kmer_size))
    {
        std::cout << kmer << ": " << code << " " << std::endl;
        assert(core::encode_kmer<core::no_ambiguity_policy>(kmer) == code);
    }
}

void encode_ambiguous_string()
{
    using namespace core;

    std::cout << "\nIteration: " << std::endl;

    auto read = std::string{ "T--ATAWAABTGTNCAAA.TT-----TTY" };
    const size_t kmer_size = 4;

    /// TODO: This is not the best way to delete gaps, since it requires changing the input
    /// string and copy data. We need a transparent gap-skipping iterator over a string/string_view
    /// that should be taken as an input argument in the core::encode/core::to_kmers functions
    read.erase(std::remove(read.begin(), read.end(), '-'), read.end());

    /// Iterate over k-mers of a string_view (implicitly created from long_read).
    /// This is more effective than calling core::encode_kmer for every k-mer of a string,
    /// because core::to_kmers calculates codes in a rolling fashion, if possible.

    /// WARNING: the amibiguous version of core::to_kmers returns a pair <ambiguous k-mer, vector of codes>
    /// It does not return resolved strings that correspond to the codes, because it requires copying
    /// of strings (there is no way to return it as std::string_view).
    size_t position = 0;
    for (const auto& [kmer, codes] : to_kmers<one_ambiguity_policy>(read, kmer_size))
    {
        std::cout << position << " " << kmer << ": ";
        for (const auto code : codes)
        {
            std::cout << code << " ";
        }
        std::cout << std::endl;

        ++position;
    }
}

int main()
{
    /// An example of core::encode_kmer and core::decode_kmer for std::string as input
    encode_unambiguous_string();

    /// An example of core::encode_kmer for std::string_view as input, and a iteration
    /// over an input seqeuence
    encode_unambiguous_string_view();

    /// An example of core::encode_kmer with one_ambiguity_policy
    encode_ambiguous_string();
}