#ifndef RAPPAS_CPP_PHYML_H
#define RAPPAS_CPP_PHYML_H

#include <string>
#include <unordered_map>
#include <memory>
#include "row.h"

namespace xpas
{
    class proba_matrix;
    class phylo_tree;


    namespace cli
    {
        class parameters;
    }

    namespace ar
    {
        /// Ancestral Reconstruction software uses their own name conventions
        /// to name internal nodes of the tree. This type maps xpas's internal node names to them
        using mapping = std::unordered_map<std::string, std::string>;

        /// Supported tools for ancestral reconstruction
        enum class software
        {
            PHYML,
            RAXML_NG
        };

        /// Ancestral reconstruction models
        enum class model
        {
            JC69,
            K80,
            F81,
            F84,
            HKY85,
            TN93,
            GTR,
            LG,
            WAG,
            JTT,
            DAYHOFF,
            DCMUT,
            CPREV,
            MTMAM,
            MTREV,
            MTART
        };

        /// All parameters for ancestral reconstruction
        struct parameters
        {
            std::string ar_dir;
            std::string binary_file;
            std::string tree_file;
            std::string alignment_file;
            model ar_model;
            double alpha;
            int categories;
        };

        /// Get necessary AR parameters from all parameters
        std::pair<ar::software, ar::parameters> make_parameters(const cli::parameters& parameters,
                                                                const std::string& ext_tree_file,
                                                                const std::string& alignment_phylip);

        /// Run ancestral reconstruction
        std::tuple<proba_matrix, phylo_tree> ancestral_reconstruction(ar::software software,
                                                                      const ar::parameters& parameters);

        /// Maps node labels of the extended tree to the node labels of the AR tree
        /// This mapping is needed to query proba_matrix.
        ar::mapping map_nodes(const phylo_tree& extended_tree, const phylo_tree& ar_tree);

    }

}



#endif