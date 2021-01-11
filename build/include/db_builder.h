#ifndef RAPPAS_CPP_DB_BUILDER_H
#define RAPPAS_CPP_DB_BUILDER_H

#include <string>
#include <xpas/phylo_kmer_db.h>
#include "extended_tree.h"
#include "ar.h"

namespace xpas
{
    class alignment;
    class proba_matrix;

    enum class filter_type
    {
        no_filter,
        entropy,
        random
    };

    xpas::phylo_kmer_db build(std::string working_directory,
                              const alignment& original_alignment, const alignment& extended_alignment,
                              const phylo_tree& original_tree, const phylo_tree& extended_tree,
                              const proba_matrix& matrix,
                              const ghost_mapping& mapping, const ar::mapping& ar_mapping,
                              bool merge_branches,
                              size_t kmer_size,
                              xpas::phylo_kmer::score_type omega,
                              filter_type filter,
                              double mu,
                              size_t num_threads);
}

#endif