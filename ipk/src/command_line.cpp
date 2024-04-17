#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "command_line.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

//--------------------------------------------------------------------------
namespace ipk::cli
{
    static std::string HELP = "help", HELP_SHORT = "h";

    /// Input files
    static std::string WORKING_DIR = "workdir", WORKING_DIR_SHORT = "w";
    static std::string OUTPUT_FILENAME = "output", OUTPUT_FILENAME_SHORT = "o";
    static std::string REFALIGN = "refalign";
    static std::string REFTREE = "reftree", REFTREE_SHORT = "t";

    /// Ancestral Reconstruction
    static std::string AR_DIR = "ar-dir";
    static std::string AR_BINARY = "ar-binary";
    static std::string AR_MODEL = "model", AR_MODEL_SHORT = "m";
    static std::string AR_ALPHA = "alpha", AR_ALPHA_SHORT = "a";
    static std::string AR_CATEGORIES = "categories";
    static std::string AR_ONLY = "ar-only";
    static std::string AR_PARAMETERS = "ar-parameters";

    /// Main options
    static std::string REDUCTION_RATIO = "reduction-ratio";
    static std::string NO_REDUCTION = "no-reduction";
    static std::string K = "k", K_SHORT = "k";
    static std::string OMEGA="omega";
    static std::string NUM_THREADS = "num-threads", NUM_THREADS_SHORT = "j";

    /// Filtering options
    static std::string MU = "mu", MU_SHORT = "u";
    static std::string MIF0 = "mif0";
    static std::string RANDOM = "random";
    static std::string MERGE_BRANCHES = "merge-branches";
    static std::string USE_UNROOTED = "use-unrooted";
    static std::string UNCOMPRESSED = "uncompressed";
    static std::string BB = "BB";
    static std::string DC = "DC";
    static std::string DCLA = "DCLA";
    static std::string DCCW = "DCCW";

    /// Ghost nodes strategy
    static std::string GHOSTS_INNER = "inner-only";
    static std::string GHOSTS_OUTER = "outer-only";
    static std::string GHOSTS_BOTH = "both";

    static std::string ON_DISK = "on-disk";


    /// Algorithm flags
    bool bb_flag = false;
    bool dc_flag = true;
    bool dcla_flag = false;
    bool dccw_flag = false;

    /// Filters flags
    bool mif0_flag = true;
    bool random_filter_flag = false;

    /// Flags for other options
    bool merge_branches_flag = false;
    bool use_unrooted_flag = false;
    bool no_reduction_flag = false;
    bool ar_only_flag = false;
    bool uncompressed_flag = false;

    /// Flags for ghost node strategy
    bool inner_only_flag = false;
    bool outer_only_flag = false;
    bool both_flag = true;

    /// Filtering algorithm flags
    bool on_disk_flag = false;

    po::options_description get_opt_description()
    {
        po::options_description desc("General options");
        desc.add_options()
            ((HELP + "," + HELP_SHORT).c_str(),
                "Show help")
            ((WORKING_DIR + "," + WORKING_DIR_SHORT).c_str(), po::value<fs::path>()->default_value(fs::current_path()),
                "Path to the working directory")
            ((OUTPUT_FILENAME + "," + OUTPUT_FILENAME_SHORT).c_str(), po::value<fs::path>()->default_value(""),
             "Output filename")
            (REFALIGN .c_str(), po::value<fs::path>()->required(),
                "Reference alignment in fasta format."
                "It must be the multiple alignment from which the reference tree was built.")
            ((REFTREE + "," + REFTREE_SHORT).c_str(), po::value<fs::path>()->required(),
                "Original phylogenetic tree file")
            (AR_DIR.c_str(), po::value<std::string>()->default_value(""),
                "Skips ancestral sequence reconstruction uses outputs from the specified directory.")
            (AR_BINARY.c_str(), po::value<std::string>()->required(),
                "Binary file for ancestral reconstruction software (PhyML, RAxML-NG).")
            ((AR_MODEL + "," + AR_MODEL_SHORT).c_str(), po::value<std::string>()->default_value("GTR"),
                "Model used in AR, one of the following:"
                "nucl  : JC69, HKY85, K80, F81, TN93, GTR"
                "amino : LG, WAG, JTT, Dayhoff, DCMut, CpREV, mMtREV, MtMam, MtArt")
            ((AR_ALPHA + "," + AR_ALPHA_SHORT).c_str(), po::value<double>()->default_value(1.0),
                "Gamma shape parameter, used in ancestral reconstruction.")
            (AR_CATEGORIES.c_str(), po::value<int>()->default_value(4),
                "Number of relative substitution rate categories, used in ancestral reconstruction.")
            ((AR_ONLY).c_str(), po::bool_switch(&ar_only_flag))
            (AR_PARAMETERS.c_str(), po::value<std::string>()->default_value(""),
                "Whitespace-separated list of arguments passed to the ancestral reconstruction tool.")

            ((K + "," + K_SHORT).c_str(), po::value<size_t>()->default_value(8),
                "k-mer length used at DB build")
            (REDUCTION_RATIO.c_str(), po::value<double>()->default_value(0.99),
                "Ratio for alignment reduction, e.g. sites holding >X% gaps are ignored.")
            (NO_REDUCTION.c_str(), po::bool_switch(&no_reduction_flag),
                "Disable alignment reduction. This will keep all sites of the reference alignment and "
                "may produce erroneous ancestral k-mers.")
            ((OMEGA).c_str(), po::value<i2l::phylo_kmer::score_type>()->default_value(1.5f),
                "Score threshold parameter")
            ((NUM_THREADS + "," + NUM_THREADS_SHORT).c_str(), po::value<size_t>()->default_value(1),
                "Number of threads")
            ((MERGE_BRANCHES).c_str(), po::bool_switch(&merge_branches_flag))
            ((USE_UNROOTED).c_str(), po::bool_switch(&use_unrooted_flag))

            ((MIF0).c_str(), po::bool_switch(&mif0_flag))
            ((RANDOM).c_str(), po::bool_switch(&random_filter_flag))
            ((MU + "," + MU_SHORT).c_str(), po::value<double>()->default_value(0.8))
            ((UNCOMPRESSED).c_str(), po::bool_switch(&uncompressed_flag))
            ((BB).c_str(), po::bool_switch(&bb_flag))
            ((DC).c_str(), po::bool_switch(&dc_flag))
            ((DCLA).c_str(), po::bool_switch(&dcla_flag))
            ((DCCW).c_str(), po::bool_switch(&dccw_flag))

            ((GHOSTS_INNER).c_str(), po::bool_switch(&inner_only_flag))
            ((GHOSTS_OUTER).c_str(), po::bool_switch(&outer_only_flag))
            ((GHOSTS_BOTH).c_str(), po::bool_switch(&both_flag))

            ((ON_DISK).c_str(), po::bool_switch(&on_disk_flag))
            ;
        return desc;
    }

    std::string get_option_list()
    {
        std::stringstream ss;
        ss << get_opt_description();
        return ss.str();
    }

    const char* to_cstr(const std::string &s)
    {
        return s.c_str();
    }

    struct cli_args
    {
        cli_args(int argc, const char* argv[])
            : argc(argc)
        {
            for (int i = 0; i < argc; ++i)
            {
                this->argv.emplace_back(argv[i]);
            }
        }

        void escape_quotes()
        {
            int start_arg = -1;
            int end_arg = -1;

            for (int i = 0; i < argc; ++i)
            {
                if ((argv[i].find('"') != std::string::npos) ||
                    (argv[i].find('\'') != std::string::npos))
                {
                    // if this is not the first quote
                    if (start_arg != -1)
                    {
                        end_arg = i;
                    }
                    else
                    {
                        start_arg = i;
                    }
                }
            }

            /// if found a match
            if (start_arg != end_arg)
            {
                argc -= end_arg - start_arg;

                /// combine args inside quotes
                std::stringstream ss;
                for (int i = start_arg; i <= end_arg; ++i)
                {
                    ss << argv[i] << " ";
                }
                argv[start_arg] = ss.str();
                argv.erase(argv.begin() + start_arg + 1, argv.begin() + end_arg + 1);
            }
        }

        [[nodiscard]]
        std::pair<int, std::vector<const char*>> to_cargs() const
        {
            std::vector<const char*>  vc;
            std::transform(argv.begin(), argv.end(), std::back_inserter(vc), to_cstr);
            return { argc, vc };
        }

        int argc;
        std::vector<std::string> argv;
    };

    parameters process_command_line(int argc, const char* argv[])
    {
        parameters parameters;
        try
        {
            /// Command line arguments can contain arguments like:
            ///     --ar-parameters "--param1 value1 --param2 value2"
            /// In this cases, quotes are ignored and parsed into many
            /// independent arguments instead of one key-value pair.
            /// We need to parse them manually
            auto args = cli_args(argc, argv);
            args.escape_quotes();
            const auto& [new_argc, new_argv] = args.to_cargs();
            /// I am sorry for this syntax
            const char* const * new_argv_array = new_argv.data();

            const po::options_description desc = get_opt_description();
            po::variables_map vm;

            po::store(po::parse_command_line(new_argc, new_argv_array, desc), vm);
            po::notify(vm);

            if (vm.count(HELP))
            {
                parameters.action = action_t::help;
                return parameters;
            }
            else
            {
                parameters.action = action_t::build;
            }

            parameters.working_directory = vm[WORKING_DIR].as<fs::path>().string();
            parameters.output_filename = vm[OUTPUT_FILENAME].as<fs::path>().string();

            /// Default output name if not given
            if (parameters.output_filename.empty())
            {
                parameters.output_filename = (parameters.working_directory / fs::path("DB.ipk")).string();
            }

            parameters.alignment_file = vm[REFALIGN].as<fs::path>().string();
            parameters.original_tree_file = vm[REFTREE].as<fs::path>().string();
            parameters.kmer_size = vm[K].as<size_t>();

            parameters.ar_dir = vm[AR_DIR].as<std::string>();
            parameters.ar_binary_file = vm[AR_BINARY].as<std::string>();
            parameters.ar_model = vm[AR_MODEL].as<std::string>();
            parameters.ar_alpha = vm[AR_ALPHA].as<double>();
            parameters.ar_categories = vm[AR_CATEGORIES].as<int>();
            parameters.ar_only = ar_only_flag;
            parameters.ar_parameters = vm[AR_PARAMETERS].as<std::string>();

            parameters.reduction_ratio = vm[REDUCTION_RATIO].as<double>();
            parameters.omega = vm[OMEGA].as<i2l::phylo_kmer::score_type>();
            parameters.num_threads = vm[NUM_THREADS].as<size_t>();
            parameters.mu = vm[MU].as<double>();

            parameters.merge_branches = merge_branches_flag;
            parameters.use_unrooted = use_unrooted_flag;
            parameters.no_reduction = no_reduction_flag;

            /// filters
            parameters.mif0_filter = mif0_flag;
            parameters.random_filter = random_filter_flag;

            /// algorithms
            parameters.bb = bb_flag;
            parameters.dc = dc_flag;
            parameters.dcla = dcla_flag;
            parameters.dccw = dccw_flag;

            parameters.uncompressed = uncompressed_flag;

            /// ghost node strategy
            parameters.inner_only = inner_only_flag;
            parameters.outer_only = outer_only_flag;
            parameters.both = both_flag;

            parameters.on_disk = on_disk_flag;
        }
        catch (const po::error& e)
        {
            throw std::runtime_error(e.what());
        }
        return parameters;
    }
}

