#ifndef RAPPAS_CPP_COMMAND_LINE_H
#define RAPPAS_CPP_COMMAND_LINE_H

#include <exception>
#include <string>
#include <map>
#include <xpas/phylo_kmer.h>
#include <pk_compute.h>

namespace xpas::cli
{
    enum class action_t
    {
        build = 0,
        help = 2
    };

    struct parameters
    {
        action_t action;

        /// Input parameters
        std::string working_directory;
        std::string alignment_file;
        std::string original_tree_file;

        // Ancestral Reconstruction parameters
        std::string ar_dir;
        std::string ar_binary_file;
        std::string ar_model;
        double ar_alpha;
        int ar_categories;
        bool ar_only;

        /// Main parameters
        double reduction_ratio;
        bool no_reduction;
        size_t kmer_size;
        xpas::phylo_kmer::score_type omega;
        size_t num_threads;

        bool merge_branches;
        bool use_unrooted;

        /// k-mer filtering parameters, mutually exclusive
        bool no_filter;
        bool entropy_filter;
        bool mif1_filter;
        bool mif0_filter;
        bool random_filter;

        /// Phylo-k-mer computation algorithms, mutually exclusize
        bool bb;
        bool dc;
        bool dcla;
        bool dccw;

        // k-mer filtering threshold
        double mu;

        // store databases uncompressed
        bool uncompressed;


        // k-mer computation algorithm
        xpas::algorithm algorithm;
    };

    std::string get_option_list();
    parameters process_command_line(int argc, const char* argv[]);
}

#endif