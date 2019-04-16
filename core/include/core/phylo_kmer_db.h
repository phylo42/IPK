#ifndef RAPPAS_CORE_PHYLO_KMER_DB_H
#define RAPPAS_CORE_PHYLO_KMER_DB_H

#include <flat_hash_map/flat_hash_map.hpp>
#include "phylo_kmer.h"

namespace core
{
    /// \brief Phylo-kmer database class, that stores all the phylo-kmers.
    class phylo_kmer_db
    {
    public:
        /// \brief A storage of a phylo-kmer information.
        /// \details Note that phylo-kmers are not stored as objects of phylo_kmer,
        /// which is just a temporary storage for a phylo-kmer information.
        using inner_storage = ska::flat_hash_map<phylo_kmer::branch_type, phylo_kmer::score_type>;
        using storage = ska::flat_hash_map<phylo_kmer::key_type, inner_storage>;
        using const_iterator = storage::const_iterator;

        phylo_kmer_db() = default;
        phylo_kmer_db(const phylo_kmer_db&) = delete;
        phylo_kmer_db(phylo_kmer_db&&) = delete;
        phylo_kmer_db& operator=(const phylo_kmer_db&) = delete;
        phylo_kmer_db& operator=(phylo_kmer_db&&) = delete;
        ~phylo_kmer_db() noexcept = default;

        /// \brief Puts a phylo-kmer in the database.
        /// \details Here we assume that all the parameters are small enough to be passed by value.
        void put(phylo_kmer::key_type key, phylo_kmer::branch_type branch, phylo_kmer::score_type score);

        const_iterator begin() const;
        const_iterator end() const;

        size_t size() const;
    private:
        storage _map;
    };
}


#endif
