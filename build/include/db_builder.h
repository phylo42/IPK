#ifndef RAPPAS_CPP_DB_BUILDER_H
#define RAPPAS_CPP_DB_BUILDER_H

#include <string>
#include <xpas/phylo_kmer_db.h>

namespace xpas
{
    class alignment;

    enum class filter_type
    {
        no_filter,
        entropy,
        random
    };

    xpas::phylo_kmer_db build(std::string working_directory,
                              std::string ar_probabilities_file,
                              xpas::alignment alignment,
                              xpas::alignment extended_alignment,
                              xpas::phylo_tree original_tree,
                              xpas::phylo_tree extended_tree,
                              bool merge_branches,
                              size_t kmer_size,
                              xpas::phylo_kmer::score_type omega,
                              filter_type filter,
                              double mu,
                              size_t num_threads);
}

#endif