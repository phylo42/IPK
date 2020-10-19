#ifndef XPAS_NEWICK_H
#define XPAS_NEWICK_H

#include <string>
#include <string_view>

namespace xpas
{
    class phylo_tree;
}

namespace xpas::io
{
    /// \brief Loads a phylogenetic tree from a newick formatted file.
    xpas::phylo_tree load_newick(const std::string& file_name);

    /// \brief Parses a phylogenetic tree from a newick formatted string.
    xpas::phylo_tree parse_newick(std::string_view newick_string);

    /// \brief Constructs a newick-formatted string from a tree (post-order traversal)
    std::string to_newick(const xpas::phylo_tree& tree);
}

#endif //RAPPAS_CORE_NEWICK_H
