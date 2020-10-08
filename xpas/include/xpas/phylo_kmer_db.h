#ifndef XPAS_PHYLO_KMER_DB_H
#define XPAS_PHYLO_KMER_DB_H


#ifdef USE_SKA_FLAT_HASH_MAP
#include <flat_hash_map/flat_hash_map.hpp>
#elif USE_SKA_BYTELL_HASH_MAP
#include <flat_hash_map/bytell_hash_map.hpp>
#elif USE_ABSL_FLAT_HASH_MAP
#include <absl/container/flat_hash_map.h>
#elif USE_FOLLY_F14_FAST_MAP
#include <folly/container/F14Map.h>
#include <folly/hash/Hash.h>
#elif USE_PHMAP_FLAT_HASH_MAP
#include <parallel_hashmap/phmap.h>
#elif USE_TSL_ROBIN_MAP
#include <tsl/robin_map.h>
#elif USE_TSL_HOPSCOTCH_MAP
#include <tsl/hopscotch_map.h>
#endif

#include "phylo_kmer.h"


namespace xpas
{
    template<typename... Args>
#ifdef USE_SKA_FLAT_HASH_MAP
    using hash_map = ska::flat_hash_map<Args...>;
#elif USE_SKA_BYTELL_HASH_MAP
    using hash_map = ska::bytell_hash_map<Args...>;
#elif USE_ABSL_FLAT_HASH_MAP
    using hash_map = absl::flat_hash_map<Args...>;
#elif USE_FOLLY_F14_FAST_MAP
    using hash_map = folly::F14FastMap<Args...>;
#elif USE_PHMAP_FLAT_HASH_MAP
    using hash_map = phmap::flat_hash_map<Args...>;
#elif USE_TSL_ROBIN_MAP
    using hash_map = tsl::robin_map<Args...>;
#elif USE_TSL_HOPSCOTCH_MAP
    using hash_map = tsl::hopscotch_map<Args...>;
#endif
    /// requires manual intervention (does not support noexcept destructor and can not be defined this way)
    ///using hash_map = robin_hood::unordered_flat_map<Args...>;
    ///using hash_map = robin_hood::unordered_map<Args...>;

    namespace impl {
        template <typename PhyloKmer>
        class search_result;

        template <>
        class search_result<unpositioned_phylo_kmer>;
    }

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

    /// \brief A phylo-kmer database that stores all phylo-kmers.
    class _base_phylo_kmer_db
    {
    public:
        /// Member types
        using key_type = phylo_kmer::key_type;


        /// Ctors, dtor and operator=
        _base_phylo_kmer_db(size_t kmer_size, xpas::phylo_kmer::score_type omega, const std::string& tree)
            : _kmer_size{ kmer_size }, _omega{ omega }, _tree { tree }
        {}

        _base_phylo_kmer_db(const _base_phylo_kmer_db&) noexcept = delete;
        _base_phylo_kmer_db(_base_phylo_kmer_db&&) = default;
        _base_phylo_kmer_db& operator=(const _base_phylo_kmer_db&) = delete;
        _base_phylo_kmer_db& operator=(_base_phylo_kmer_db&&) = default;
        virtual ~_base_phylo_kmer_db() noexcept = default;

        /// \brief Returns the k-mer size.
        [[nodiscard]]
        size_t kmer_size() const noexcept
        {
            return _kmer_size;
        }

        void set_kmer_size(size_t kmer_size) noexcept
        {
            _kmer_size = kmer_size;
        }

        /// \brief Returns omega (the parameter of xpas::score_threshold)
        [[nodiscard]]
        phylo_kmer::score_type omega() const noexcept
        {
            return _omega;
        }

        void set_omega(phylo_kmer::score_type omega)
        {
            _omega = omega;
        }

        /// \brief Returns a view to the newick formatted phylogenetic tree
        [[nodiscard]]
        std::string_view tree() const noexcept
        {
            return _tree;
        }

        void set_tree(const std::string& tree)
        {
            _tree = tree;
        }

    private:

        /// \brief K-mer size.
        /// \details This number is given by user to the constructor. We can not guarantee
        /// that the keys stored in hash tables actually correspond to substrings of the size _kmer_size.
        /// Example: DNA ('A', 'C', 'G', 'T")
        ///     key('AAA') == key('AA') == 0
        /// e.g. putting 0 in hashtable, we assume it corresponds to 'AAA' having _kmer_size == 3,
        /// but we can not guarantee that it was not calculated for another k-mer size by mistake.
        size_t _kmer_size;

        /// \brief Score threshold paramenter.
        /// \sa core::score_threshold
        xpas::phylo_kmer::score_type _omega;

        /// \brief Newick formatted phylogenetic tree
        std::string _tree;
    };

    namespace impl
    {
        /// \brief A search result wrapper around a collection of pairs [branch, score]
        /// to iterate over.

        template <typename ValueType>
        class search_result
        {
        public:
            using const_iterator = typename ValueType::const_iterator;

            search_result() noexcept = default;
            search_result(const_iterator begin, const_iterator end) noexcept
                : _begin{ begin }, _end{ end }
            {}

            search_result(const search_result&) noexcept = default;
            ~search_result() noexcept = default;

            [[nodiscard]]
            const_iterator begin() const noexcept
            {
                return _begin;
            }

            [[nodiscard]]
            const_iterator end() const noexcept
            {
                return _end;
            }
        private:
            const_iterator _begin;
            const_iterator _end;
        };
    }

    template <typename PhyloKmer>
    class _phylo_kmer_db : public _base_phylo_kmer_db
    {};
}

namespace xpas {
    template <>
    class _phylo_kmer_db<unpositioned_phylo_kmer> : public _base_phylo_kmer_db
    {
    public:
        using value_type = std::vector<_pkdb_value<false>>;

        /// \brief A storage of a phylo-kmer information.
        /// \details Note that phylo-kmers are not stored as objects of phylo_kmer,
        /// which is just a temporary storage for a phylo-kmer information.
        using storage = hash_map<key_type, value_type>;
        using const_iterator = typename storage::const_iterator;

        _phylo_kmer_db(size_t kmer_size, xpas::phylo_kmer::score_type omega, const std::string& tree)
        : _base_phylo_kmer_db(kmer_size, omega, tree)
        {}
        _phylo_kmer_db(const _phylo_kmer_db&) noexcept = delete;
        _phylo_kmer_db(_phylo_kmer_db&&) = default;
        _phylo_kmer_db& operator=(const _phylo_kmer_db&) = delete;
        _phylo_kmer_db& operator=(_phylo_kmer_db&&) = default;
        ~_phylo_kmer_db() noexcept override = default;

        /// Access
        /// \brief Searches for a key against the database.
        /// \details WARNING: This method does not know how the key was calculated. It is required
        /// to provide keys of substrings of size _kmer_size to get correct results.
        /// \sa _kmer_size
        [[nodiscard]]
        std::optional<impl::search_result<value_type>> search(key_type key) const noexcept
        {
            if (auto it = _map.find(key); it != _map.end())
            {
                return impl::search_result<value_type>{ it->second.begin(), it->second.end() };
            }
            else
            {
                return std::nullopt;
            }
        }

        /// Iterators
        /// \brief Returns an iterator to the beginning
        [[nodiscard]]
        const_iterator begin() const noexcept
        {
            return std::begin(_map);
        }

        /// \brief Returns an iterator to the end
        [[nodiscard]]
        const_iterator end() const noexcept
        {
            return std::end(_map);
        }

        /// Capacity
        /// \brief Returns the number of keys
        [[nodiscard]]
        size_t size() const noexcept
        {
            return _map.size();
        }

        /// \brief Returns a hash function used to hash kmers
        [[nodiscard]]
        typename storage::hasher hash_function() const noexcept
        {
            return _map.hash_function();
        }

        /// Modifiers
        /// \brief Puts a phylo-kmer information in the database.
        /// \details This method is unsafe, which means it does not control if the value has
        /// a correct branch id, the score is maximal etc. All of this checks must be done before calling
        /// this method. It just puts the value in a hash map.
        /// WARNING: This method does not know how the key was calculated. Here we assume it represents
        /// a string of size _kmer_size.
        /// \sa _kmer_size
        void unsafe_insert(key_type key, const unpositioned_pkdb_value& value)
        {
            _map[key].push_back(value);
        }

    private:
        storage _map;
    };

    template <>
    class _phylo_kmer_db<positioned_phylo_kmer> : public _base_phylo_kmer_db
    {
    public:
        using pkdb_value_type = positioned_pkdb_value;
        using value_type = std::vector<positioned_pkdb_value>;

        /// \brief A storage of a phylo-kmer information.
        /// \details Note that phylo-kmers are not stored as objects of phylo_kmer,
        /// which is just a temporary storage for a phylo-kmer information.
        using storage = hash_map<key_type, value_type>;
        using const_iterator = typename storage::const_iterator;

        _phylo_kmer_db(size_t kmer_size, xpas::phylo_kmer::score_type omega, const std::string& tree)
            : _base_phylo_kmer_db(kmer_size, omega, tree)
        {}
        _phylo_kmer_db(const _phylo_kmer_db&) noexcept = delete;
        _phylo_kmer_db(_phylo_kmer_db&&) = default;
        _phylo_kmer_db& operator=(const _phylo_kmer_db&) = delete;
        _phylo_kmer_db& operator=(_phylo_kmer_db&&) = default;
        ~_phylo_kmer_db() noexcept override = default;

        /// Access
        /// \brief Searches for a key against the database.
        /// \details WARNING: This method does not know how the key was calculated. It is required
        /// to provide keys of substrings of size _kmer_size to get correct results.
        /// \sa _kmer_size
        [[nodiscard]]
        std::optional<impl::search_result<value_type>> search(key_type key) const noexcept
        {
            if (auto it = _map.find(key); it != _map.end())
            {
                return impl::search_result<value_type>{ it->second.begin(), it->second.end() };
            }
            else
            {
                return std::nullopt;
            }
        }

        /// Iterators
        /// \brief Returns an iterator to the beginning
        [[nodiscard]]
        const_iterator begin() const noexcept
        {
            return std::begin(_map);
        }

        /// \brief Returns an iterator to the end
        [[nodiscard]]
        const_iterator end() const noexcept
        {
            return std::end(_map);
        }

        /// Capacity
        /// \brief Returns the number of keys
        [[nodiscard]]
        size_t size() const noexcept
        {
            return _map.size();
        }

        /// \brief Returns a hash function used to hash kmers
        [[nodiscard]]
        typename storage::hasher hash_function() const noexcept
        {
            return _map.hash_function();
        }

        /// Modifiers
        /// \brief Puts a phylo-kmer information in the database.
        /// \details This method is unsafe, which means it does not control if the value has
        /// a correct branch id, the score is maximal etc. All of this checks must be done before calling
        /// this method. It just puts the value in a hash map.
        /// WARNING: This method does not know how the key was calculated. Here we assume it represents
        /// a string of size _kmer_size.
        /// \sa _kmer_size
        void unsafe_insert(key_type key, const pkdb_value_type& value)
        {
            _map[key].push_back(value);
        }

        /// \brief Replace a phylo-kmer in the database.
        /// \details This method will remove all values associated to the given key,
        /// and will add a new value.
        void replace(key_type key, const pkdb_value_type& value)
        {
            _map[key].clear();
            _map[key].push_back(value);
        }

    private:
        storage _map;
    };


    using phylo_kmer_db = _phylo_kmer_db<phylo_kmer>;
}

#endif
