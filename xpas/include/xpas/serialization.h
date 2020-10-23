#ifndef XPAS_SERIALIZATION_H
#define XPAS_SERIALIZATION_H

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <fstream>
#include "phylo_kmer_db.h"
#include "version.h"

namespace xpas
{
    /// A struct with definitions for different version of the serialization protocol
    struct protocol
    {
        /// xpas v0.1.x
        static const int v0_1_x = 1;

        /// xpas v0.2.x
        static const int v0_2_WITHOUT_POSITIONS = 2;
        static const int v0_2_WITH_POSITIONS = 3;

#ifdef KEEP_POSITIONS
        static const unsigned int CURRENT = v0_2_WITH_POSITIONS;
#else
        static const unsigned int CURRENT = v0_2_WITHOUT_POSITIONS;
#endif
    };

    xpas::phylo_kmer_db load(const std::string& filename)
    {
        std::ifstream ifs(filename);
        ::boost::archive::binary_iarchive ia(ifs);

        xpas::phylo_kmer_db db { 0, 0.0, "", "" };
        ia & db;
        return db;
    }

    void save(const xpas::phylo_kmer_db& db, const std::string& filename)
    {
        std::ofstream ofs(filename);
        ::boost::archive::binary_oarchive oa(ofs);
        oa & db;
    }
}

namespace boost::serialization
{
    template<class Archive>
    inline void save(Archive& ar,
                     const xpas::_phylo_kmer_db<xpas::positioned_phylo_kmer>& db,
                     unsigned int /*version*/)
    {
        ar & std::string(db.sequence_type());

        const auto original_tree_view = std::string{ db.tree() };
        ar & original_tree_view;

        const auto kmer_size = db.kmer_size();
        ar & kmer_size;

        xpas::phylo_kmer::score_type omega = db.omega();
        ar & omega;

        size_t table_size = db.size();
        ar & table_size;

        for (const auto& [key, entries] : db)
        {
            size_t entries_size = entries.size();
            ar & key & entries_size;

            for (const auto& [branch, score, position] : entries)
            {
                ar & branch & score & position;
            }
        }
    }

    template<class Archive>
    inline void save(Archive& ar,
                     const xpas::_phylo_kmer_db<xpas::unpositioned_phylo_kmer>& db,
                     unsigned int /*version*/)
    {
        ar & std::string(db.sequence_type());

        const auto original_tree_view = std::string{ db.tree() };
        ar & original_tree_view;

        size_t kmer_size = db.kmer_size();
        ar & kmer_size;

        xpas::phylo_kmer::score_type omega = db.omega();
        ar & omega;

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
    inline void load(Archive& ar,
                     xpas::_phylo_kmer_db<xpas::positioned_phylo_kmer>& db,
                     const unsigned int version)
    {
        /// Early versions are not supported
        if (version < xpas::protocol::v0_2_WITH_POSITIONS)
        {
            throw std::runtime_error("Failed to load database: the database does not have positional information.");
        }

        /// if we did not throw an exception, positions will be loaded.
        db.set_positions_loaded(true);

        /// Deserialization of the content added in versions v0.2.x
        if (version > xpas::protocol::v0_1_x)
        {
            std::string sequence_type;
            ar & sequence_type;
            db.set_sequence_type(std::move(sequence_type));
        }

        /// Deserialization of the main content
        {
            std::string tree;
            ar & tree;
            db.set_tree(std::move(tree));

            size_t kmer_size = 0;
            ar & kmer_size;
            db.set_kmer_size(kmer_size);

            xpas::phylo_kmer::score_type omega = 0;
            ar & omega;
            db.set_omega(omega);

            size_t table_size = 0;
            ar & table_size;
            for (size_t i = 0; i < table_size; ++i)
            {
                xpas::phylo_kmer::key_type key = xpas::phylo_kmer::na_key;
                size_t entries_size = 0;
                ar & key;
                ar & entries_size;
                for (size_t j = 0; j < entries_size; ++j)
                {
                    xpas::phylo_kmer::branch_type branch = xpas::phylo_kmer::na_branch;
                    xpas::phylo_kmer::score_type score = xpas::phylo_kmer::na_score;
                    xpas::phylo_kmer::pos_type position = xpas::phylo_kmer::na_pos;
                    ar & branch & score & position;
                    db.unsafe_insert(key, { branch, score, position });
                }
            }
        }
    }

    template<class Archive>
    inline void load(Archive& ar,
                     xpas::_phylo_kmer_db<xpas::unpositioned_phylo_kmer>& db,
                     const unsigned int version)
    {
        /// Early versions are not supported
        if (version < xpas::protocol::v0_1_x)
        {
            throw std::runtime_error("Failed to load database: this database was built with older version of xpas.");
        }

        db.set_positions_loaded(false);

        /// Deserialization of the content added in versions v0.2.x
        if (version > xpas::protocol::v0_1_x)
        {
            std::string sequence_type;
            ar & sequence_type;
            db.set_sequence_type(std::move(sequence_type));
        }

        std::string tree;
        ar & tree;
        db.set_tree(std::move(tree));

        size_t kmer_size = 0;
        ar & kmer_size;
        db.set_kmer_size(kmer_size);

        xpas::phylo_kmer::score_type omega = 0;
        ar & omega;
        db.set_omega(omega);

        size_t table_size = 0;
        ar & table_size;
        for (size_t i = 0; i < table_size; ++i)
        {
            xpas::phylo_kmer::key_type key = xpas::phylo_kmer::na_key;
            size_t entries_size = 0;
            ar & key;
            ar & entries_size;
            for (size_t j = 0; j < entries_size; ++j)
            {
                /// classic deserialization of non-positioned phylo k-mers
                xpas::phylo_kmer::branch_type branch = xpas::phylo_kmer::na_branch;
                xpas::phylo_kmer::score_type score = xpas::phylo_kmer::na_score;
                ar & branch & score;

                /// if the database has positions, read and ignore them
                if (version == xpas::protocol::v0_2_WITH_POSITIONS)
                {
                    xpas::phylo_kmer::pos_type position = xpas::phylo_kmer::na_pos;
                    ar & position;
                }
                db.unsafe_insert(key, { branch, score });
            }
        }
    }

    // split non-intrusive serialization function member into separate
    // non intrusive save/load member functions
    template<class Archive>
    inline void serialize(Archive& ar, xpas::phylo_kmer_db& db, unsigned int file_version)
    {
        boost::serialization::split_free(ar, db, file_version);
    }
}

BOOST_CLASS_VERSION(xpas::phylo_kmer_db, xpas::protocol::CURRENT)

#endif //RAPPAS_CORE_SERIALIZATION_H
