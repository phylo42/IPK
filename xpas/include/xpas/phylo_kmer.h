#ifndef XPAS_PHYLO_KMER_H
#define XPAS_PHYLO_KMER_H

#include "seq.h"
#include <limits>
#include <string>
#include <optional>

namespace xpas
{
    /// \brief A phylo k-mer structure.
    /// \details A key-value pair for a phylo-kmer, where key is a key_type value of a k-mer, and value is
    /// posterior probability score of this k-mer. Branch node id, position etc. omitted here, because these values
    /// are shared among multiple phylo k-mers and can be stored more effectively.
    //template <bool KeepPositions> struct _phylo_kmer
    //    : std::conditional<KeepPositions, impl::_positioned_phylo_kmer, impl::_unpositioned_phylo_kmer>::type
    struct unpositioned_phylo_kmer
    {
        /// The same type used to store a k-mer value
        /// (essentially we do not distinguish between "k-mer value" and "phylo-kmer value")
        using key_type = seq_traits::key_type;

        /// The type of a "posterior probability" score of a phylokmer
        using score_type = float;

        /// The type of a branch node id, that a phylokmer is mapped to.
        using branch_type = uint32_t;

        /// The type of a phylokmer's position in the alignment
        using pos_type = uint16_t;

        static constexpr key_type na_key = std::numeric_limits<key_type>::max();
        static constexpr score_type na_score = std::numeric_limits<score_type>::quiet_NaN();
        static constexpr branch_type na_branch = std::numeric_limits<branch_type>::max();
        static constexpr pos_type na_pos = std::numeric_limits<pos_type>::max();

        [[nodiscard]]
        bool is_nan() const
        {
            return (key == na_key) && (score == na_score);
        }

        key_type key;
        score_type score;
    };

    struct positioned_phylo_kmer : public unpositioned_phylo_kmer
    {
        positioned_phylo_kmer(key_type key, score_type score, pos_type pos)
            : unpositioned_phylo_kmer{ key, score }, position{ pos }
        {}

        pos_type position;
    };


#ifdef KEEP_POSITIONS
    using phylo_kmer = positioned_phylo_kmer;
#else
    using phylo_kmer = unpositioned_phylo_kmer;
#endif

    /// \brief A value of phylo-kmer stored in phylo_kmer_db.
    /// \details One k-mer can correspond to multiple values, so the key of a k-mer is not stored here.
    template <bool KeepPositions>
    struct _pkdb_value;

    template <>
    struct _pkdb_value<true>
    {
        phylo_kmer::branch_type branch;
        phylo_kmer::score_type score;
        phylo_kmer::pos_type position;

        _pkdb_value(phylo_kmer::branch_type _branch, phylo_kmer::score_type _score,
                    phylo_kmer::pos_type _position)
            : branch{ _branch }
            , score{ _score }
            , position{ _position }
        {}
    };

    template <>
    struct _pkdb_value<false>
    {
        phylo_kmer::branch_type branch;
        phylo_kmer::score_type score;

        _pkdb_value(phylo_kmer::branch_type _branch, phylo_kmer::score_type _score)
            : branch{ _branch }, score{ _score } {}

    };

    using positioned_pkdb_value = _pkdb_value<true>;
    using unpositioned_pkdb_value = _pkdb_value<false>;

#ifdef KEEP_POSITIONS
    using pkdb_value = positioned_pkdb_value;
#else
    using pkdb_value = unpositioned_pkdb_value;
#endif

    /// A simple struct to connect (un)positioned versions of phylo_kmer and pkdb_value types.
    template<typename PhyloKmer>
    struct get_pkdb_value_type
    {
    };

    template <>
    struct get_pkdb_value_type<unpositioned_phylo_kmer>
    {
        using type = unpositioned_pkdb_value;
    };

    template <>
    struct get_pkdb_value_type<positioned_phylo_kmer>
    {
        using type = positioned_pkdb_value;
    };

    /// A templated factory function that unifies the interface of creating pkdb_values with
    /// and without positions.
    template <typename PhyloKmerType,
              typename PkDbValue = typename get_pkdb_value_type<PhyloKmerType>::type>
    PkDbValue make_pkdb_value(phylo_kmer::branch_type branch, phylo_kmer::score_type score,
                              phylo_kmer::pos_type position);

    /// A templated factory function that unifies the interface of creating phylo k-mers with
    /// and without positions.
    template <typename PhyloKmerType>
    PhyloKmerType make_phylo_kmer(phylo_kmer::key_type key, phylo_kmer::score_type score,
                                  phylo_kmer::pos_type position);


    bool operator==(const positioned_phylo_kmer& lhs, const positioned_phylo_kmer& rhs) noexcept;
    bool operator==(const unpositioned_phylo_kmer& lhs, const unpositioned_phylo_kmer& rhs) noexcept;

    /// Returns a minumum score
    phylo_kmer::score_type score_threshold(phylo_kmer::score_type omega, size_t kmer_size);

    /// \brief Returns one or more codes of the input k-mer (depending on the policy)
    /// \details Assumes that the size of input sequence equals k
    template<typename AmbiguityPolicy>
    std::optional<typename AmbiguityPolicy::value_type> encode_kmer(std::string_view kmer);

    template<typename AmbiguityPolicy>
    std::optional<phylo_kmer::key_type> encode_kmer(const std::string& kmer)
    {
        return encode_kmer<AmbiguityPolicy>(std::string_view{ kmer });
    }

    /// Creates a string of size kmer_size by given key
    std::string decode_kmer(phylo_kmer::key_type key, size_t kmer_size);
}

#endif