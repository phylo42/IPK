#include <iostream>
#include <sstream>
#include <chrono>
#include <boost/filesystem.hpp>
#include <xpas/phylo_kmer_db.h>
#include <xpas/serialization.h>
#include "phylo_tree.h"
#include "command_line.h"
#include "exceptions.h"
#include "db_builder.h"
#include "alignment.h"
#include "return.h"
#include "ar.h"
#include "proba_matrix.h"

namespace fs = boost::filesystem;

return_code print_help()
{
    std::cout << "RAPPAS2" << std::endl << std::endl
              << "Usage: rappas2 [...]" << std::endl
              << cli::get_option_list() << std::endl;

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

void check_parameters(const cli::cli_parameters& parameters)
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
    xpas::save_alignment(alignment, fasta_path.string());

    fs::path phylip_path = directory / "extended_align.phylip";
    //xpas::save_alignment(alignment, phylip_path.string());


    return { fasta_path.string(), phylip_path.string() };
}

return_code build_database(const cli::cli_parameters& parameters)
{
    if (parameters.kmer_size > xpas::seq_traits::max_kmer_length)
    {
        std::cerr << "Maximum k-mer size allowed: " << xpas::seq_traits::max_kmer_length << std::endl;
        return return_code::argument_error;
    }

    std::cout << "TODO: MEASURE TIME HERE" << std::endl << std::endl;
    auto alignment = xpas::preprocess_alignment(parameters.working_directory,
                                                parameters.alignment_file,
                                                parameters.reduction_ratio,
                                                parameters.no_reduction);

    auto [original_tree, extended_tree, mapping] = xpas::preprocess_tree(parameters.original_tree_file);
    const auto extended_tree_file = save_extended_tree(parameters.working_directory, extended_tree);

    auto extended_alignment = xpas::extend_alignment(alignment, extended_tree);
    const auto& [ext_alignment_fasta, ext_alignment_phylip] = save_extended_alignment(parameters.working_directory, extended_alignment);

    const auto proba_matrix = xpas::ancestral_reconstruction(extended_tree_file, ext_alignment_phylip);

    const auto db = xpas::build(parameters.working_directory,
                                parameters.ar_probabilities_file,
                                std::move(alignment),
                                std::move(extended_alignment),
                                std::move(original_tree),
                                std::move(extended_tree),
                                parameters.merge_branches, parameters.kmer_size, parameters.omega,
                                xpas::filter_type::no_filter, parameters.mu,
                                parameters.num_threads);

    const auto db_filename = fs::path(parameters.working_directory) / generate_db_name(db);

    std::cout << "Saving database to: " << db_filename.string() << "..." << std::endl;
    const auto begin = std::chrono::steady_clock::now();
    xpas::save(db, db_filename.string());
    std::cout << "Time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - begin).count() << std::endl << std::endl;

    return return_code::success;
}

return_code run(const cli::cli_parameters& parameters)
{
    switch (parameters.action)
    {
        case cli::help:
        {
            return print_help();
        }
        case cli::build:
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
        const cli::cli_parameters parameters = cli::process_command_line(argc, argv);
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