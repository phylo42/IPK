#ifndef RAPPAS_CORE_PHYLO_KMER_DB_H
#define RAPPAS_CORE_PHYLO_KMER_DB_H

#include <flat_hash_map/flat_hash_map.hpp>
#include "phylo_kmer.h"

namespace core
{
    namespace impl {
        class search_result;
    }

    /// \brief Phylo-kmer database class, that stores all the phylo-kmers.
    class phylo_kmer_db
    {
    public:
        using key_type = phylo_kmer::key_type;
        using inner_key_type = phylo_kmer::branch_type;
        using value_type = phylo_kmer::score_type;

        /// \brief A storage of a phylo-kmer information.
        /// \details Note that phylo-kmers are not stored as objects of phylo_kmer,
        /// which is just a temporary storage for a phylo-kmer information.
        using inner_storage = ska::flat_hash_map<inner_key_type, value_type>;
        using storage = ska::flat_hash_map<key_type, inner_storage>;
        using const_iterator = storage::const_iterator;

        phylo_kmer_db() = default;
        phylo_kmer_db(const phylo_kmer_db&) = delete;
        phylo_kmer_db(phylo_kmer_db&&) noexcept = default;
        phylo_kmer_db& operator=(const phylo_kmer_db&) = delete;
        phylo_kmer_db& operator=(phylo_kmer_db&&) noexcept = default;
        ~phylo_kmer_db() noexcept = default;

        /// \brief Puts a phylo-kmer in the database.
        /// \details Here we assume that all the parameters are small enough to be passed by value.
        void put(key_type key, inner_key_type branch, value_type score);

        /// \brief Searches for a key against the database.
        std::optional<impl::search_result> search(key_type key) const;

        const_iterator begin() const;
        const_iterator end() const;

        size_t size() const;
    private:
        storage _map;
    };

    namespace impl
    {
        /// \brief A wrapper around a collection of pairs [branch, score] for a search result
        /// to iterate over.
        class search_result
        {
        public:
            using const_iterator = phylo_kmer_db::inner_storage::const_iterator;

            search_result() noexcept;
            search_result(const_iterator begin, const_iterator end) noexcept;
            search_result(const search_result&) noexcept = default;
            ~search_result() noexcept = default;

            const_iterator begin() const;
            const_iterator end() const;

        private:
            const_iterator _begin;
            const_iterator _end;
        };
    }
}


#endif
