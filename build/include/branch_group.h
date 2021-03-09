#ifndef XPAS_BRANCH_GROUP_H
#define XPAS_BRANCH_GROUP_H

#include <xpas/hash_map.h>
#include <xpas/phylo_kmer.h>
#include <xpas/phylo_kmer_db.h>

namespace xpas
{

    /// \brief A hash map to store all the phylo-kmers, placed to one original node
#ifdef KEEP_POSITIONS

    struct score_pos_pair
        {
            phylo_kmer::score_type score;
            phylo_kmer::pos_type position;
        };

        using group_hash_map = hash_map<phylo_kmer::key_type, score_pos_pair>;
#else
    using group_hash_map = hash_map<phylo_kmer::key_type, phylo_kmer::score_type>;
#endif

    void save_group_map(const group_hash_map& map, const std::string& filename);

    group_hash_map load_group_map(const std::string& filename);

    std::string get_groups_dir(const std::string& working_dir);

    /// \brief Returns a filename of the group hashmap
    std::string get_group_map_file(const std::string& working_dir,
                                   const phylo_kmer::branch_type& group, size_t batch_idx);

    /// Merges hashmaps of the same index into a database
    xpas::phylo_kmer_db merge_batch(const std::string& working_dir,
                                    const std::vector<phylo_kmer::branch_type>& group_ids, size_t batch_idx);

    /// Puts a kmer in the hash. Takes a maximum score between the existing value
    /// of the k-mer (if any) and the provided value.
    void put(group_hash_map& map, const phylo_kmer& kmer);

    void accumulate(group_hash_map& map, const phylo_kmer& kmer,
                    phylo_kmer::score_type default_value, phylo_kmer::score_type log_rev_threshold);

    /// Returns the number of the batch of the given k-mer
    size_t kmer_batch(phylo_kmer::key_type key, size_t n_ranges);

}

#endif //XPAS_BRANCH_GROUP_H
