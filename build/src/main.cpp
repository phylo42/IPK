#include <iostream>
#include <sstream>
#include <chrono>
#include <boost/filesystem.hpp>
#include <xpas/phylo_kmer_db.h>
#include <xpas/serialization.h>
#include <xpas/phylo_tree.h>
#include "command_line.h"
#include "exceptions.h"
#include "db_builder.h"
#include "extended_tree.h"
#include "alignment.h"
#include "return.h"
#include "ar.h"
#include "proba_matrix.h"
#include "filter.h"

namespace fs = boost::filesystem;

return_code print_help()
{
    std::cout << "RAPPAS2" << std::endl << std::endl
              << "Usage: rappas2 [...]" << std::endl
              << xpas::cli::get_option_list() << std::endl;

    return return_code::help;
}

std::string generate_db_name(const xpas::phylo_kmer_db& db)
{
    const auto kmer_size = db.kmer_size();
    const auto omega = db.omega();

    /// RAPPAS has to support the following output name convention for omega:
    /// o0.5,
    /// o0.75
    /// o1.25
    /// o1.12345
    /// o1.5
    ///
    /// BUT 2.0 not 2
    ///     1.0 not 1
    ///
    /// Do not ask me why...

    /// std::setprecision can not handle this case.
    /// So we have to a trailing zero value for round values of omega
    std::string omega_str = std::to_string(omega);
    const auto last_zero = omega_str.find_last_not_of('0') + 1;
    omega_str.erase(last_zero, omega_str.size() - last_zero);
    if (static_cast<int>(omega) == omega) {
        omega_str += "0";
    }

    std::ostringstream out;
    out << "DB_k" << kmer_size << "_o" << omega_str << ".rps";
    return out.str();
}

void check_parameters(const xpas::cli::parameters& parameters)
{
   if (!xpas::keep_positions && parameters.merge_branches)
   {
       throw std::runtime_error("--merge-branches is only supported for xpas compiled with the KEEP_POSITIONS flag.");
   }
}

std::string save_extended_tree(const std::string& working_dir, const xpas::phylo_tree& tree)
{
    fs::path directory = fs::path(working_dir) / "extended_trees";
    fs::path full_path = directory / "extended_tree.newick";
    fs::create_directories(directory);
    xpas::save_tree(tree, full_path.string());
    return full_path.string();
}

std::pair<std::string, std::string> save_extended_alignment(const std::string& working_dir, const xpas::alignment& alignment)
{
    fs::path directory = fs::path(working_dir) / "extended_trees";
    fs::create_directories(directory);


    fs::path fasta_path = directory / "extended_align.fasta";
    std::cout << "Saving alignment to " << fasta_path.string() << "..." << std::endl;
    xpas::save_alignment(alignment, fasta_path.string(), xpas::alignment_format::FASTA);

    fs::path phylip_path = directory / "extended_align.phylip";
    std::cout << "Saving alignment to " << phylip_path.string() << "..." << std::endl;
    xpas::save_alignment(alignment, phylip_path.string(), xpas::alignment_format::PHYLIP);

    return { fasta_path.string(), phylip_path.string() };
}

std::string save_rerooted_tree(const std::string& working_dir, const xpas::phylo_tree& tree)
{
    fs::path directory = fs::path(working_dir) / "AR";
    fs::create_directories(directory);

    fs::path tree_path = directory / "ar_tree_rerooted.newick";
    std::cout << "Saving tree to " << tree_path.string() << "..." << std::endl;
    xpas::save_tree(tree, tree_path.string());
    return tree_path.string();
}

xpas::filter_type get_filter_type(const xpas::cli::parameters& parameters)
{
    if (parameters.entropy_filter)
    {
        return xpas::filter_type::entropy;
    }
    else if (parameters.mis_filter)
    {
        return xpas::filter_type::mis;
    }
    else if (parameters.mif_filter)
    {
        return xpas::filter_type::mif;
    }
    else if (parameters.random_filter)
    {
        return xpas::filter_type::random;
    }

    return xpas::filter_type::no_filter;
}

return_code build_database(const xpas::cli::parameters& parameters)
{
    if (parameters.kmer_size > xpas::seq_traits::max_kmer_length)
    {
        std::cerr << "Maximum k-mer size allowed: " << xpas::seq_traits::max_kmer_length << std::endl;
        return return_code::argument_error;
    }

    /// Load and filter the reference alignment
    auto alignment = xpas::preprocess_alignment(parameters.working_directory,
                                                parameters.alignment_file,
                                                parameters.reduction_ratio,
                                                parameters.no_reduction);

    /// Load and extend the reference tree
    const auto& [original_tree, extended_tree, ghost_mapping] = xpas::preprocess_tree(parameters.original_tree_file,
                                                                                      parameters.use_unrooted);
    const auto extended_tree_file = save_extended_tree(parameters.working_directory, extended_tree);

    /// Extend the alignment
    auto extended_alignment = xpas::extend_alignment(alignment, extended_tree);
    const auto& [ext_alignment_fasta, ext_alignment_phylip] =
        save_extended_alignment(parameters.working_directory, extended_alignment);

    /// Prepare and run ancestral reconstruction
    auto [ar_software, ar_parameters] = xpas::ar::make_parameters(parameters,
                                                                  extended_tree_file, ext_alignment_phylip);
    auto [proba_matrix, ar_tree] = xpas::ar::ancestral_reconstruction(ar_software, ar_parameters);

    if (parameters.ar_only)
    {
        std::cout << "--ar-only requested. Finishing after ancestral reconstruction." << std::endl;
        return return_code::success;
    }

    /// Ancestral Reconstruction will unroot the input tree even if it is rooted:
    /// ((A, B)node, C)root  ===>   (C, A, B)newick_root;
    if (original_tree.is_rooted() && !ar_tree.is_rooted())
    {
        /// Re-root it back:
        /// ((A,B)newick_root,C)added_root;
        xpas::reroot_tree(ar_tree);
        save_rerooted_tree(parameters.working_directory, ar_tree);
    }

    const auto ar_mapping = xpas::ar::map_nodes(extended_tree, ar_tree);

    /// Generate phylo k-mers
    const auto db = xpas::build(parameters.working_directory,
                                original_tree, extended_tree,
                                proba_matrix,
                                ghost_mapping, ar_mapping,
                                parameters.merge_branches, parameters.kmer_size, parameters.omega,
                                get_filter_type(parameters), parameters.mu,
                                parameters.num_threads);

    /// Deserialize database
    const auto db_filename = fs::path(parameters.working_directory) / generate_db_name(db);
    std::cout << "Saving database to: " << db_filename.string() << "..." << std::endl;
    const auto begin = std::chrono::steady_clock::now();
    xpas::save(db, db_filename.string());
    std::cout << "Time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - begin).count() << std::endl << std::endl;

    return return_code::success;
}

return_code run(const xpas::cli::parameters& parameters)
{
    switch (parameters.action)
    {
        case xpas::cli::action_t::help:
        {
            return print_help();
        }
        case xpas::cli::action_t::build:
        {
            return build_database(parameters);
        }
        default:
        {
            return return_code::unknown_error;
        }
    }
}

int main(int argc, const char* argv[])
{
    std::ios::sync_with_stdio(false);

    try
    {
        const auto parameters = xpas::cli::process_command_line(argc, argv);
        check_parameters(parameters);
        run(parameters);
    }
    catch (const conflicting_options& e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (const bad_options& e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}