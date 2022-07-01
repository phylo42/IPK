#include "branch_group.h"
#include <iostream>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/filesystem.hpp>

using namespace xpas;
namespace fs = boost::filesystem;

void xpas::save_group_map(const group_hash_map& map, const std::string& filename)
{
    std::ofstream ofs(filename);
    boost::archive::binary_oarchive oa(ofs);
    oa & map;
}

/// \brief Loads a hash map from file
group_hash_map xpas::load_group_map(const std::string& filename)
{
    std::ifstream ifs(filename);
    boost::archive::binary_iarchive ia(ifs);

    group_hash_map map;
    ia & map;
    return map;
}

std::string xpas::get_groups_dir(const std::string& working_dir)
{
    return { (fs::path{working_dir} / fs::path{"hashmaps"}).string() };
}

std::string xpas::get_group_map_file(const std::string& working_dir, const phylo_kmer::branch_type& group, size_t batch_idx)
{
    return (get_groups_dir(working_dir) /
            fs::path{std::to_string(group) + "_" + std::to_string(batch_idx) + ".hash"}).string();
}

phylo_kmer_db xpas::merge_batch(const std::string& working_dir,
                                const std::vector<phylo_kmer::branch_type>& group_ids, size_t batch_idx)
{
    //std::cout << "Merging hash maps [batch index = " << batch_idx << "]..." << std::endl;
    phylo_kmer_db temp_db(0, 1.0, xpas::seq_type::name, "");

    /// Load hash maps and merge them
    for (const auto group_id : group_ids)
    {
        const auto hash_map = load_group_map(get_group_map_file(working_dir, group_id, batch_idx));
#ifdef KEEP_POSITIONS
        for (const auto& [key, score_pos_pair] : hash_map)
            {
                const auto& [score, position] = score_pos_pair;
                temp_db.unsafe_insert(key, {group_id, score, position});
            }
#else
        for (const auto& [key, score] : hash_map)
        {
            temp_db.unsafe_insert(key, {group_id, score});
        }
#endif
    }

    return temp_db;
}

#ifdef KEEP_POSITIONS
void xpas::put(group_hash_map& map, const phylo_kmer& kmer)
    {
        if (auto it = map.find(kmer.key); it != map.end())
        {
            if (it->second.score < kmer.score)
            {
                map[kmer.key] = { kmer.score, kmer.position };
            }
        }
        else
        {
            map[kmer.key] = { kmer.score, kmer.position };
        }
    }
#else
void xpas::put(group_hash_map& map, const phylo_kmer& kmer)
{
    if (auto it = map.find(kmer.key); it != map.end())
    {
        if (it->second < kmer.score)
        {
            map[kmer.key] = kmer.score;
        }
    }
    else
    {
        map[kmer.key] = kmer.score;
    }
}
#endif

size_t xpas::kmer_batch(phylo_kmer::key_type key, size_t n_ranges)
{
    return key % n_ranges;
}

namespace boost::serialization
{

#ifdef KEEP_POSITIONS
    /// Serialize a hash map
    template<class Archive>
    inline void save(Archive& ar, const group_hash_map& map, const unsigned int /*version*/)
    {
        size_t map_size = map.size();
        ar & map_size;

        for (const auto& [key, score_pos_pair] : map)
        {
            const auto& [score, position] = score_pos_pair;
            ar & key & score & position;
        }
    }

    /// Deserialize a hash map
    template<class Archive>
    inline void load(Archive& ar, group_hash_map& map, unsigned int version)
    {
        (void)version;

        size_t map_size = 0;
        ar & map_size;

        for (size_t i = 0; i < map_size; ++i)
        {
            xpas::phylo_kmer::key_type key = xpas::phylo_kmer::na_key;
            xpas::phylo_kmer::score_type score = xpas::phylo_kmer::na_score;
            xpas::phylo_kmer::pos_type position = xpas::phylo_kmer::na_pos;
            ar & key & score & position;
            map[key] = { score, position };
        }
    }
#else
    /// Serialize a hash map
    template<class Archive>
    inline void save(Archive& ar, const group_hash_map& map, const unsigned int /*version*/)
    {
        size_t map_size = map.size();
        ar & map_size;

        for (const auto&[key, score] : map)
        {
            ar & key & score;
        }
    }

    /// Deserialize a hash map
    template<class Archive>
    inline void load(Archive& ar, group_hash_map& map, unsigned int version)
    {
        (void)version;

        size_t map_size = 0;
        ar & map_size;

        for (size_t i = 0; i < map_size; ++i)
        {
            xpas::phylo_kmer::key_type key = xpas::phylo_kmer::na_key;
            xpas::phylo_kmer::score_type score = xpas::phylo_kmer::na_score;
            ar & key & score;
            map[key] = score;
        }
    }
#endif

    // split non-intrusive serialization function member into separate
    // non intrusive save/load member functions
    template<class Archive>
    inline void serialize(Archive& ar, group_hash_map& map, unsigned int file_version)
    {
        boost::serialization::split_free(ar, map, file_version);
    }
}
