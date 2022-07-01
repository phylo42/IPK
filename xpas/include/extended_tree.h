#ifndef XPAS_TREE_H
#define XPAS_TREE_H

#include <xcl/phylo_kmer.h>
#include <xcl/phylo_tree.h>


namespace xpas
{
    using ghost_mapping = std::unordered_map<std::string, xpas::phylo_kmer::branch_type>;

    /// Read and preprocess a phylogentic tree
    std::tuple<phylo_tree, phylo_tree, ghost_mapping> preprocess_tree(const std::string& filename, bool use_unrooted);

    /// Reroot tree if needed
    /// changes the tree from (a, b, c); to ((b, c), a);
    void reroot_tree(phylo_tree& tree);
}

#endif //XPAS_TREE_H
