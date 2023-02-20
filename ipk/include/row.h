#ifndef RAPPAS_CPP_ROW_H
#define RAPPAS_CPP_ROW_H

#include <i2l/phylo_kmer.h>
#include <array>

namespace ipk
{
    using i2l::phylo_kmer;

    struct proba_pair
    {
        phylo_kmer::score_type score;
        phylo_kmer::key_type index;
    };

    using branch_type = phylo_kmer::branch_type;
    using row_type = std::array<proba_pair, i2l::seq_traits::alphabet_size>;
}

#endif //RAPPAS_CPP_ROW_H
