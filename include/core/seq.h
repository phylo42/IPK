#ifndef RAPPAS_CORE_SEQ_H
#define RAPPAS_CORE_SEQ_H

#include <cstdint>
#include <cstddef>

namespace core
{
    /// \brief Auxiliary structure to define the seq_type.
    /// \sa seq_type
    struct dna
    {};
    struct protein
    {};

    /// \brief Sequence traits type.
    /// \details Describes all the properties based on the size of an alphabet, the alphabet itself etc.
    ///
    /// Must declare:
    ///     - size_t alphabet_size
    ///     - size_t max_kmer_length
    ///     - char_type decode(size_t)
    ///     - size_t encode(char_type)
    template<typename SeqType>
    struct seq_traits_impl;

    template<>
    struct seq_traits_impl<dna>
    {
        /// \brief The type used to store one base of a sequence.
        using char_type = uint8_t;

        /// \brief The type used to store a k-mer value.
        /// \details K-mers are not stored as strings or char* values, but as values of this type instead.
        /// For example for DNA and k==3 the k-mer "AAA" == 000ul if kmer_t is unsigned long.
        using key_type = uint64_t;

        static constexpr char_type char_set[] = {'A', 'C', 'G', 'T'};

        static constexpr char_type decode(size_t code)
        {
            return char_set[code];
        }

        static constexpr size_t encode(char_type base)
        {
            switch (base)
            {
                case 'A':
                    return 0;
                case 'C':
                    return 1;
                case 'G':
                    return 2;
                case 'T':
                    return 3;
                default:
                    return 4;
            }
        }

        static constexpr char_type ambiguous_chars[] =  {'N', '.', '-', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'};
        static constexpr size_t alphabet_size = sizeof(char_set);

        static constexpr size_t max_kmer_length = 16;
    };

    /// \brief A compile-time constant for a current sequence type.
    /// \details Sequence type is determined at compile time for efficiency reasons. To use another sequence type
    /// (e.g. proteins), declare a new seq_traits_impl and recompile core for a new sequence type.
    using seq_traits = seq_traits_impl<dna>;

    /// \brief Returns amount of bits used to store one base of a sequence of given type.
    template<typename SeqType>
    constexpr seq_traits::key_type bit_length();

    template<>
    constexpr seq_traits::key_type bit_length<dna>()
    {
        return 2ul;
    }

    /// \brief Returns a value of kmer_t type that can mask the rightmost base of a kmer_t.
    /// \details E.g. for DNA 0b1111...1100
    template<typename SeqType>
    constexpr seq_traits::key_type rightest_symbol_mask();

    template<>
    constexpr seq_traits::key_type rightest_symbol_mask<dna>()
    {
        return ~0b11ul;
    }
}

#endif