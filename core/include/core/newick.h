#ifndef RAPPAS_CORE_NEWICK_H
#define RAPPAS_CORE_NEWICK_H

#include <stack>
#include <string>

namespace core
{
    class phylo_tree;
}

namespace rappas
{
    namespace io
    {
        /// \brief Loads a phylogenetic tree from a newick formatted file.
        core::phylo_tree load_newick(const std::string& file_name);
    }
}

#endif //RAPPAS_CORE_NEWICK_H
