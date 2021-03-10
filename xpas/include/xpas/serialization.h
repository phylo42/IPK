#ifndef XPAS_SERIALIZATION_H
#define XPAS_SERIALIZATION_H

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include "phylo_kmer_db.h"
#include "version.h"

namespace fs = boost::filesystem;

namespace xpas
{
    /// A struct with definitions for different version of the serialization protocol
    struct protocol
    {
        /// xpas v0.1.x
        static const int v0_1_x = 2;

        /// xpas v0.2.x
        static const int v0_2_WITHOUT_POSITIONS = 3;
        static const int v0_2_WITH_POSITIONS = 4;

        /// The protocol is the same for v0.2.x - v0.3.0
        static const int v0_3_0 = 5;

        /// xpas v0.3.1+
        static const int v0_3_1_WITHOUT_POSITIONS = 5;
        static const int v0_3_1_WITH_POSITIONS = 5;


#ifdef KEEP_POSITIONS
        static const unsigned int CURRENT = v0_3_1_WITH_POSITIONS;
#else
        static const unsigned int CURRENT = v0_3_1_WITHOUT_POSITIONS;
#endif
    };

    xpas::phylo_kmer_db load_compressed(const std::string& filename)
    {
        std::ifstream ifs(filename);
        boost::iostreams::filtering_istream in;
        boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
        in.push(boost::iostreams::zlib_decompressor(zp));
        in.push(ifs);

        boost::archive::binary_iarchive ia(in);

        xpas::phylo_kmer_db db { 0, 0.0, "", "", score_model_type::MAX};
        ia & db;
        return db;
    }

    xpas::phylo_kmer_db load_uncompressed(const std::string& filename)
    {
        std::ifstream ifs(filename);
        boost::archive::binary_iarchive ia(ifs);

        xpas::phylo_kmer_db db { 0, 0.0, "", "", score_model_type::MAX};
        ia & db;
        return db;
    }


    xpas::phylo_kmer_db load(const std::string& filename)
    {
        if (!fs::exists(filename))
        {
            throw std::runtime_error("No such file: " + filename);
        }

        /// Versions earlier than v0.2.1 were not using zlib compression.
        /// There is no good way to figure out if the file is compressed or not
        /// than to just try to decompress first.
        try
        {
            return load_compressed(filename);
        }
        catch (const boost::iostreams::zlib_error& error)
        {
            return load_uncompressed(filename);
        }
    }

    void save(const xpas::phylo_kmer_db& db, const std::string& filename)
    {
        std::ofstream ofs(filename);

        boost::iostreams::filtering_ostream out;
        boost::iostreams::zlib_params zp(boost::iostreams::zlib::best_speed);
        out.push(boost::iostreams::zlib_compressor(zp));
        out.push(ofs);

        ::boost::archive::binary_oarchive oa(out);
        oa & db;
    }
}

namespace boost::serialization
{
    template<class Archive>
    inline void save(Archive& ar,
                     const xpas::_phylo_kmer_db<xpas::positioned_phylo_kmer>& db,
                     unsigned int version)
    {
        ar & std::string(db.sequence_type());

        const auto original_tree_view = std::string{ db.tree() };
        ar & original_tree_view;

        const auto kmer_size = db.kmer_size();
        ar & kmer_size;

        xpas::phylo_kmer::score_type omega = db.omega();
        ar & omega;

        if (version > xpas::protocol::v0_3_0)
        {
            xpas::score_model_type score_model = db.score_model();
            ar & score_model;

            xpas::phylo_kmer::score_type threshold = db.threshold();
            ar & threshold;
        }

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
                     unsigned int version)
    {
        ar & std::string(db.sequence_type());

        const auto original_tree_view = std::string{ db.tree() };
        ar & original_tree_view;

        size_t kmer_size = db.kmer_size();
        ar & kmer_size;

        xpas::phylo_kmer::score_type omega = db.omega();
        ar & omega;

        if (version > xpas::protocol::v0_3_0)
        {
            xpas::score_model_type score_model = db.score_model();
            ar & score_model;

            xpas::phylo_kmer::score_type threshold = db.threshold();
            ar & threshold;
        }

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

            /// Deserialization of the content added in v0.3.1+
            if (version > xpas::protocol::v0_3_0)
            {
                xpas::score_model_type score_model = xpas::score_model_type::MAX;
                ar & score_model;
                db.set_score_model(score_model);

                xpas::phylo_kmer::score_type threshold = 0.0f;
                ar & threshold;
                db.set_threshold(threshold);
            }

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

        /// Deserialization of the content added in v0.3.1+
        if (version > xpas::protocol::v0_3_0)
        {
            xpas::score_model_type score_model = xpas::score_model_type::MAX;
            ar & score_model;
            db.set_score_model(score_model);

            xpas::phylo_kmer::score_type threshold = 0.0f;
            ar & threshold;
            db.set_threshold(threshold);
        }

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
