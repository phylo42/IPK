#ifndef IPK_DB_BUILDER_H
#define IPK_DB_BUILDER_H

#include <string>
#include <i2l/phylo_kmer_db.h>
#include "extended_tree.h"
#include "ar.h"

namespace i2l
{
    class phylo_tree;
}

namespace ipk
{
    class alignment;
    class proba_matrix;
    enum class filter_type;
    enum class algorithm;
    enum class ghost_strategy;

    void build(const std::string& working_directory, const std::string& output_filename,
               const i2l::phylo_tree& original_tree, const i2l::phylo_tree& extended_tree,
               proba_matrix& matrix,
               const ghost_mapping& mapping, const ar::mapping& ar_mapping, bool merge_branches,
               ipk::algorithm algorithm, ipk::ghost_strategy strategy,
               size_t kmer_size, i2l::phylo_kmer::score_type omega,
               filter_type filter, double mu,
               size_t num_threads, bool on_disk);
}

#endif