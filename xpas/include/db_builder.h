#ifndef RAPPAS_CPP_DB_BUILDER_H
#define RAPPAS_CPP_DB_BUILDER_H

#include <string>
#include <xcl/phylo_kmer_db.h>
#include "extended_tree.h"
#include "ar.h"

namespace xcl
{
    class phylo_tree;
}

namespace xpas
{
    class alignment;
    class proba_matrix;
    enum class filter_type;
    enum class algorithm;

    xcl::phylo_kmer_db build(const std::string& working_directory,
                             const xcl::phylo_tree& original_tree, const xcl::phylo_tree& extended_tree,
                             const proba_matrix& matrix,
                             const ghost_mapping& mapping, const ar::mapping& ar_mapping, bool merge_branches,
                             xpas::algorithm algorithm, size_t kmer_size, xpas::phylo_kmer::score_type omega,
                             filter_type filter, double mu,
                             size_t num_threads);
}

#endif