#ifndef XPAS_NEWICK_H
#define XPAS_NEWICK_H

#include <string>
#include <string_view>

namespace xcl
{
    class phylo_tree;
}

namespace xcl::io
{
    /// \brief Loads a phylogenetic tree from a newick formatted file.
    xcl::phylo_tree load_newick(const std::string& file_name);

    /// \brief Parses a phylogenetic tree from a newick formatted string.
    xcl::phylo_tree parse_newick(std::string_view newick_string);

    /// \brief Constructs a Newick-formatted string from the input tree.
    /// Depending on the jplace parameter, it builds
    ///     false) Pure newick: (label:branch_length,label:branch_length)...
    ///     true) Jplace: (label:branch_length{node_postorder_id}, ...)...
    std::string to_newick(const xcl::phylo_tree& tree, bool jplace=false);
}

/// \brief Outputs a tree in Jplace-extended Newick format:
/// (label:branch_length{node_postorder_id}, ...)
std::ostream& operator<<(std::ostream& out, const xcl::phylo_tree& tree);

#endif //RAPPAS_CORE_NEWICK_H
