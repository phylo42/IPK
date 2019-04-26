#ifndef RAPPAS_CORE_SERIALIZATION_H
#define RAPPAS_CORE_SERIALIZATION_H

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <fstream>
#include "phylo_kmer_db.h"

namespace core
{
    ::core::phylo_kmer_db load(const std::string& filename)
    {
        std::ifstream ifs(filename);
        boost::archive::binary_iarchive ia(ifs);

        size_t kmer_size;
        ia & kmer_size;
        ::core::phylo_kmer_db db { kmer_size };
        ia & db;
        return db;
    }

    void save(const ::core::phylo_kmer_db& db, const std::string& filename)
    {
        std::ofstream ofs(filename);
        boost::archive::binary_oarchive oa(ofs);
        oa & db;
    }
}

namespace boost {
    namespace serialization
    {
        template<class Archive>
        inline void save(Archive& ar, const ::core::phylo_kmer_db& db, const unsigned int /* file_version */)
        {
            size_t kmer_size = db.kmer_size();
            ar & kmer_size;

            size_t table_size = db.size();
            ar & table_size;

            for (const auto& [key, entries] : db)
            {
                size_t entries_size = entries.size();
                ar & key & entries_size;
                for (const auto& [branch, score] : entries)
                {
                    ar & branch & score;
                }
            }
        }

        template<class Archive>
        inline void load(Archive& ar, ::core::phylo_kmer_db& db, const unsigned int /* file_version */)
        {
            /// Note:
            /// k-mer size is already loaded

            size_t table_size = 0;
            ar & table_size;
            for (size_t i = 0; i < table_size; ++i)
            {
                ::core::phylo_kmer::key_type key = ::core::phylo_kmer::nan_key;
                size_t entries_size = 0;
                ar & key;
                ar & entries_size;
                for (size_t j = 0; j < entries_size; ++j)
                {
                    ::core::phylo_kmer::branch_type branch = ::core::phylo_kmer::nan_branch;
                    ::core::phylo_kmer::score_type score = ::core::phylo_kmer::nan_score;
                    ar & branch & score;
                    db.put(key, branch, score);
                }
            }
        }

        // split non-intrusive serialization function member into separate
        // non intrusive save/load member functions
        template<class Archive>
        inline void serialize(Archive& ar, ::core::phylo_kmer_db& db, const unsigned int file_version)
        {
            boost::serialization::split_free(ar, db, file_version);
        }
    }
}

#endif //RAPPAS_CORE_SERIALIZATION_H
