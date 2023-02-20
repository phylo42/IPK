#ifndef XPAS_TREE_H
#define XPAS_TREE_H

#include <i2l/phylo_kmer.h>
#include <i2l/phylo_tree.h>


namespace ipk
{
    using ghost_mapping = std::unordered_map<std::string, i2l::phylo_kmer::branch_type>;

    /// Read and preprocess a phylogentic tree
    std::tuple<i2l::phylo_tree, i2l::phylo_tree, ghost_mapping> preprocess_tree(const std::string& filename, bool use_unrooted);

    /// Reroot tree if needed
    /// changes the tree from (a, b, c); to ((b, c), a);
    void reroot_tree(i2l::phylo_tree& tree);
}

#endif //XPAS_TREE_H
