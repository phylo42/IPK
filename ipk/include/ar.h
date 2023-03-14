#ifndef AR_H
#define AR_H

#include <string>
#include <unordered_map>
#include <memory>
#include "row.h"

namespace i2l
{
    class phylo_tree;
}

namespace ipk
{
    class proba_matrix;

    namespace cli
    {
        class parameters;
    }

    namespace ar
    {
        /// Ancestral Reconstruction software uses their own name conventions
        /// to name internal nodes of the tree. This type maps IPK's internal node names to them
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
            JC,
            K80,
            F81,
            HKY,
            F84,
            TN93ef,
            TN93,
            K81,
            K81uf,
            TPM2,
            TPM2uf,
            TPM3,
            TPM3uf,
            TIM1,
            TIM1uf,
            TIM2,
            TIM2uf,
            TIM3,
            TIM3uf,
            TVMef,
            TVM,
            SYM,
            GTR,
            // FIXME: add other protein models
            JTT,
            LG,
            WAG
        };

        /// All parameters for ancestral reconstruction
        struct parameters
        {
            std::string ar_dir;
            std::string binary_file;
            std::string tree_file;
            std::string alignment_file;

            /// User-defined arguments for the AR software
            std::string ar_parameters;

            /// If ar_parameters are empty, those three are used
            model ar_model;
            double alpha;
            int categories;
        };

        /// Get necessary AR parameters from all parameters
        std::pair<ar::software, ar::parameters> make_parameters(const cli::parameters& parameters,
                                                                const std::string& ext_tree_file,
                                                                const std::string& alignment_phylip);

        /// Run ancestral reconstruction
        std::tuple<proba_matrix, i2l::phylo_tree> ancestral_reconstruction(ar::software software,
                                                                           const ar::parameters& parameters);

        /// Maps node labels of the extended tree to the node labels of the AR tree
        /// This mapping is needed to query proba_matrix.
        ar::mapping map_nodes(const i2l::phylo_tree& extended_tree, const i2l::phylo_tree& ar_tree);

    }

}



#endif