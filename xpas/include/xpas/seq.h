#ifndef XPAS_SEQ_H
#define XPAS_SEQ_H

#include <cstdint>
#include <cstddef>
#include <vector>
#include "optional.h"

namespace xpas
{
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

#ifdef SEQ_TYPE_DNA

    /// \brief Auxiliary structure to define the seq_type.
    /// \sa seq_type
    struct dna
    {
        inline static const char* name = "DNA";
    };

    /// \brief A compile-time constant for a current sequence type.
    /// \details Sequence type is determined at compile time for efficiency reasons.
    using seq_type = dna;

    template<>
    struct seq_traits_impl<dna>
    {
        /// \brief The type used to store one base of a sequence.
        using char_type = unsigned char;

        /// \brief The type used to store a k-mer value.
        /// \details K-mers are not stored as strings or char* values, but as values of this type instead.
        /// For example for DNA and k==3 the k-mer "AAA" == 000ul if kmer_t is unsigned long.
        using key_type = uint32_t;

        static const char_type code_to_key[];

        /// \brief Alphabet size
        /// \details The number of different *codes* of the alphabet. For DNA, 'T' and 'U' have the
        /// same code, which counts only once.
        static constexpr size_t alphabet_size = 4;
        static constexpr size_t max_kmer_length = 12;

        /// A type returned by encoding an unambiguous character
        using unambiguous_code_t = optional<uint8_t>;

        /// A type returned by encoding an ambiguous character
        using ambiguous_code_t = optional<std::vector<uint8_t>>;

        static constexpr unambiguous_code_t key_to_code(char_type base)
        {
            switch (base)
            {
                case 'A':
                    [[fallthrough]];
                case 'a':
                    return { 0 };
                case 'C':
                    [[fallthrough]];
                case 'c':
                    return { 1 };
                case 'G':
                    [[fallthrough]];
                case 'g':
                    return { 2 };
                case 'T':
                    [[fallthrough]];
                case 't':
                    [[fallthrough]];
                case 'U':
                    [[fallthrough]];
                case 'u':
                    return { 3 };
                case 'N':
                    [[fallthrough]];
                case 'n':
                    [[fallthrough]];
                case '-':
                    [[fallthrough]];
                case '.':
                    [[fallthrough]];
                default:
                    return nullopt;
            }
        };

        static ambiguous_code_t ambiguous_key_to_code(char_type base)
        {
            switch (base)
            {
                /// Unambiguous characters
                case 'A':
                    [[fallthrough]];
                case 'a':
                    return { { 0 } };
                case 'C':
                    [[fallthrough]];
                case 'c':
                    return { { 1 } };
                case 'G':
                    [[fallthrough]];
                case 'g':
                    return { { 2 } };
                case 'T':
                    [[fallthrough]];
                case 't':
                    [[fallthrough]];
                case 'U':
                    [[fallthrough]];
                case 'u':
                    return { { 3 } };
                /// Purine
                /// R = A | G
                case 'R':
                    [[fallthrough]];
                case 'r':
                    return { { 0, 2 } };
                /// Pyrimidine
                /// Y = C | T
                case 'Y':
                    [[fallthrough]];
                case 'y':
                    return { { 1, 3 } };
                /// Strong
                /// S = C | G
                case 'S':
                    [[fallthrough]];
                case 's':
                    return { { 1, 2 } };
                /// Weak
                /// W = A | T
                case 'W':
                    [[fallthrough]];
                case 'w':
                    return { { 0, 3 } };
                /// Keto
                /// K = G | T
                case 'K':
                    [[fallthrough]];
                case 'k':
                    return { { 2, 3 } };
                /// Amino
                /// M = A | C
                case 'M':
                    [[fallthrough]];
                case 'm':
                    return { { 0, 1 } };
                /// Not A
                case 'B':
                    [[fallthrough]];
                case 'b':
                    return { { 1, 2, 3 } };
                /// Not T/U
                case 'V':
                    [[fallthrough]];
                case 'v':
                    return { { 0, 1, 2 } };
                /// Not G
                case 'H':
                    [[fallthrough]];
                case 'h':
                    return { { 0, 1, 3 } };
                /// Not C
                case 'D':
                    [[fallthrough]];
                case 'd':
                    return { { 0, 2, 3 } };
                /// Any
                case 'X':
                    [[fallthrough]];
                case 'x':
                    [[fallthrough]];
                case 'N':
                    [[fallthrough]];
                case 'n':
                    [[fallthrough]];
                case '.':
                    return { { 0, 1, 2, 3 } };
                default:
                    return nullopt;
            }
        }

        static bool is_gap(char_type base)
        {
            return base == '-';
        }

        static bool is_ambiguous(char_type base)
        {
            const auto optional = ambiguous_key_to_code(base);
            return !optional || optional->size() > 1;
        }
    };

#elif SEQ_TYPE_AA

    /// \brief Auxiliary structure to define the seq_type.
    /// \sa seq_type
    struct aa
    {
        inline static const char* name = "Proteins";
    };

    /// \brief A compile-time constant for a current sequence type.
    /// \details Sequence type is determined at compile time for efficiency reasons.
    using seq_type = aa;

    template<>
    struct seq_traits_impl<aa>
    {
        /// \brief The type used to store one base of a sequence.
        using char_type = unsigned char;

        /// \brief The type used to store a k-mer value.
        /// \details K-mers are not stored as strings or char* values, but as values of this type instead.
        /// For example for proteins and k==3 the k-mer "RRR" == 000ul if kmer_t is unsigned long.
        using key_type = uint32_t;

        static const char_type code_to_key[];

        /// \brief Alphabet size
        static constexpr size_t alphabet_size = 20;
        static constexpr size_t max_kmer_length = 6;

        /// A type returned by encoding an unambiguous character
        using unambiguous_code_t = std::optional<uint8_t>;

        /// A type returned by encoding an ambiguous character
        using ambiguous_code_t = std::optional<std::vector<uint8_t>>;

        static constexpr unambiguous_code_t key_to_code(char_type base)
        {
            switch (base)
            {
                case 'R':
                    [[fallthrough]];
                case 'r':
                    return { 0 };
                case 'H':
                    [[fallthrough]];
                case 'h':
                    return { 1 };
                case 'K':
                    [[fallthrough]];
                case 'k':
                    return { 2 };
                case 'D':
                    [[fallthrough]];
                case 'd':
                    return { 3 };
                case 'E':
                    [[fallthrough]];
                case 'e':
                    return { 4 };
                case 'S':
                    [[fallthrough]];
                case 's':
                    return { 5 };
                case 'T':
                    [[fallthrough]];
                case 't':
                    return { 6 };
                case 'N':
                    [[fallthrough]];
                case 'n':
                    return { 7 };
                case 'Q':
                    [[fallthrough]];
                case 'q':
                    return { 8 };
                case 'C':
                    [[fallthrough]];
                case 'c':
                    [[fallthrough]];
                case 'U':
                    [[fallthrough]];
                case 'u':
                    return { 9 };
                case 'G':
                    [[fallthrough]];
                case 'g':
                    return { 10 };
                case 'P':
                    [[fallthrough]];
                case 'p':
                    return { 11 };
                case 'A':
                    [[fallthrough]];
                case 'a':
                    return { 12 };
                case 'I':
                    [[fallthrough]];
                case 'i':
                    return { 13 };
                case 'L':
                    [[fallthrough]];
                case 'l':
                    [[fallthrough]];
                case 'O':
                    [[fallthrough]];
                case 'o':
                    return { 14 };
                case 'M':
                    [[fallthrough]];
                case 'm':
                    return { 15 };
                case 'F':
                    [[fallthrough]];
                case 'f':
                    return { 16 };
                case 'W':
                    [[fallthrough]];
                case 'w':
                    return { 17 };
                case 'Y':
                    [[fallthrough]];
                case 'y':
                    return { 18 };
                case 'V':
                    [[fallthrough]];
                case 'v':
                    return { 19 };
                /// the same as in RAPPAS
                case '*':
                    [[fallthrough]];
                case '!':
                    [[fallthrough]];
                case '-':
                    [[fallthrough]];
                case '.':
                    [[fallthrough]];
                case 'X':
                    [[fallthrough]];
                case 'x':
                    [[fallthrough]];
                default:
                    return std::nullopt;
            }
        };

        static ambiguous_code_t ambiguous_key_to_code(char_type base)
        {
            switch (base)
            {
                /// Unambiguous characters
                case 'R':
                    [[fallthrough]];
                case 'r':
                    return {{0}};
                case 'H':
                    [[fallthrough]];
                case 'h':
                    return {{1}};
                case 'K':
                    [[fallthrough]];
                case 'k':
                    return {{2}};
                case 'D':
                    [[fallthrough]];
                case 'd':
                    return {{3}};
                case 'E':
                    [[fallthrough]];
                case 'e':
                    return {{4}};
                case 'S':
                    [[fallthrough]];
                case 's':
                    return {{5}};
                case 'T':
                    [[fallthrough]];
                case 't':
                    return {{6}};
                case 'N':
                    [[fallthrough]];
                case 'n':
                    return {{7}};
                case 'Q':
                    [[fallthrough]];
                case 'q':
                    return {{8}};
                case 'C':
                    [[fallthrough]];
                case 'c':
                    [[fallthrough]];
                case 'U':
                    [[fallthrough]];
                case 'u':
                    return {{9}};
                case 'G':
                    [[fallthrough]];
                case 'g':
                    return {{10}};
                case 'P':
                    [[fallthrough]];
                case 'p':
                    return {{11}};
                case 'A':
                    [[fallthrough]];
                case 'a':
                    return {{12}};
                case 'I':
                    [[fallthrough]];
                case 'i':
                    return {{13}};
                case 'L':
                    [[fallthrough]];
                case 'l':
                    [[fallthrough]];
                case 'O':
                    [[fallthrough]];
                case 'o':
                    return {{14}};
                case 'M':
                    [[fallthrough]];
                case 'm':
                    return {{15}};
                case 'F':
                    [[fallthrough]];
                case 'f':
                    return {{16}};
                case 'W':
                    [[fallthrough]];
                case 'w':
                    return {{17}};
                case 'Y':
                    [[fallthrough]];
                case 'y':
                    return {{18}};
                case 'V':
                    [[fallthrough]];
                case 'v':
                    return {{19}};

                /// Ambiguous characters
                /// B = D | N
                case 'B':
                    [[fallthrough]];
                case 'b':
                    return {{3, 7}};
                /// Z = E | Q
                case 'Z':
                    [[fallthrough]];
                case 'z':
                    return {{4, 8}};
                /// J = I | L
                case 'J':
                    [[fallthrough]];
                case 'j':
                    return {{13, 14}};

                /// Any
                case 'X':
                    [[fallthrough]];
                case 'x':
                    [[fallthrough]];
                case '*':
                    [[fallthrough]];
                case '!':
                    [[fallthrough]];
                case '.':
                    return {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}};
                default:
                    return std::nullopt;
            }
        }

        static bool is_gap(char_type base)
        {
            return base == '-';
        }

        static bool is_ambiguous(char_type base)
        {
            const auto optional = ambiguous_key_to_code(base);
            return !optional || optional->size() > 1;
        }
    };
#else
    static_assert(false, """Please define a sequence type to compile core. Supported types:\n"""
                         """SEQ_TYPE_DNA"""
                         """SEQ_TYPE_AA""");
#endif

    using seq_traits = seq_traits_impl<seq_type>;

    /// \brief Returns amount of bits used to store one base of a sequence of given type.
    template<typename SeqType>
    constexpr seq_traits::key_type bit_length();

    /// \brief Returns a value of kmer_t type that can mask the rightmost base of a kmer_t.
    /// \details E.g. for DNA 0b1111...1100
    ///               for AA 0b11...1100000
    template<typename SeqType>
    constexpr seq_traits::key_type rightmost_symbol_mask();

#ifdef SEQ_TYPE_DNA
    template<>
    constexpr seq_traits::key_type bit_length<dna>()
    {
        return seq_traits::key_type{ 2u };
    }

    template<>
    constexpr seq_traits::key_type rightmost_symbol_mask<dna>()
    {
        return seq_traits::key_type{ ~0b11u };
    }

#elif SEQ_TYPE_AA
    template<>
    constexpr seq_traits::key_type bit_length<aa>()
    {
        return seq_traits::key_type{ 5u };
    }

    template<>
    constexpr seq_traits::key_type rightmost_symbol_mask<aa>()
    {
        return seq_traits::key_type{ ~0b11111u };
    }
#else
    /// ...
#endif

    template<typename SeqTraits>
    constexpr typename SeqTraits::char_type decode_impl(uint8_t code)
    {
        return SeqTraits::code_to_key[code];
    }

    const auto decode = decode_impl<seq_traits>;

    template<class SeqTraits>
    struct no_ambiguity_policy_impl
    {
        /// a k-mer code
        using value_type = typename SeqTraits::key_type;

        /// an encode function
        constexpr static auto encode = SeqTraits::key_to_code;
    };

    template<class SeqTraits>
    struct one_ambiguity_policy_impl
    {
        /// a vector of k-mer codes
        using value_type = std::vector<typename SeqTraits::key_type>;

        /// an encode function
        constexpr static auto encode = SeqTraits::ambiguous_key_to_code;
    };

    /// Encodes a character of a sequence of current used sequence type
    template<typename AmbiguityPolicy>
    auto encode(seq_traits::char_type key)
    {
        return AmbiguityPolicy::encode(key);
    }

    /// Aliases of ambiguity polices for used sequence type
    using no_ambiguity_policy = no_ambiguity_policy_impl<seq_traits>;
    using one_ambiguity_policy = one_ambiguity_policy_impl<seq_traits>;
}

#endif