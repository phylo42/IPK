#include <iostream>
#include <numeric>
#include <cmath>
#include <string>
#include <optional>
#include <regex>
#include <array>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/process.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <csv-parser/csv.h>
#include <i2l/newick.h>
#include <i2l/phylo_tree.h>
#include "ar.h"
#include "row.h"
#include "proba_matrix.h"
#include "command_line.h"

using std::string;
using std::cout, std::endl;
namespace bp = boost::process;
namespace fs = boost::filesystem;

namespace ipk::ar
{
    /// \brief Reads ancestral reconstruction output
    class ar_reader
    {
    public:
        virtual ~ar_reader() noexcept = default;
        virtual ipk::proba_matrix read() = 0;
    };

    /// \brief Reads a PhyML output into a matrix.
    class phyml_reader : public ar_reader
    {
    public:
        phyml_reader(const string& file_name) noexcept;
        phyml_reader(const phyml_reader&) = delete;
        phyml_reader(phyml_reader&&) = delete;
        phyml_reader& operator=(const phyml_reader&) = delete;
        phyml_reader& operator=(phyml_reader&&) = delete;
        ~phyml_reader() noexcept override = default;

        proba_matrix read() override;

    private:
        proba_matrix read_matrix();

        string _file_name;
    };

    /// \brief Reads RAXML-NG output into a matrix.
    class raxmlng_reader : public ar_reader
    {
    public:
        raxmlng_reader(const string& file_name) noexcept;
        raxmlng_reader(const phyml_reader&) = delete;
        raxmlng_reader(raxmlng_reader&&) = delete;
        raxmlng_reader& operator=(const raxmlng_reader&) = delete;
        raxmlng_reader& operator=(raxmlng_reader&&) = delete;
        ~raxmlng_reader() noexcept override = default;

        proba_matrix read();

    private:
        proba_matrix read_matrix();

        string _file_name;
    };

    phyml_reader::phyml_reader(const string& file_name) noexcept
        : _file_name{ file_name }
    {}

    proba_matrix phyml_reader::read()
    {
        try
        {
            cout << "Loading PhyML results: " + _file_name << "..." << endl;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            auto matrix = read_matrix();

            cout << "Loaded " << matrix.num_branches() << " matrices of " <<
                 matrix.num_sites() << " rows." << endl;
            cout << "Time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - begin).count() << endl << endl;
            return matrix;
        }
        catch (::io::error::integer_overflow& error)
        {
            throw std::runtime_error("PhyML result parsing error: " + string(error.what()));
        }
    }

    proba_matrix phyml_reader::read_matrix()
    {
#ifdef SEQ_TYPE_DNA
        std::ios::sync_with_stdio(false);

        proba_matrix result;
        std::ifstream infile(_file_name);

        bool is_header = true;

        string line;
        while (std::getline(infile, line))
        {
            if (is_header)
            {
                // if the line starts with 'Site', the header is finished
                if (line.size() > 4 && line.rfind("Site", 0) == 0)
                {
                    is_header = false;
                }
            }
            else
            {
                size_t site = 0;
                std::string node_label;
                ipk::phylo_kmer::score_type a, c, g, t;

                std::istringstream iss(line);
                if (!(iss >> site >> node_label >> a >> c >> g >> t))
                {
                    throw std::runtime_error("Parsing error: could not parse the line " + line);
                }

                auto new_column = std::vector<phylo_kmer::score_type>{ a, c, g, t };

                /// log-transform the probabilities
                auto log = [](auto value) { return std::log10(value); };
                std::transform(begin(new_column), end(new_column), begin(new_column), log);

                auto& node_matrix = result[node_label];
                node_matrix.set_label(node_label);
                node_matrix.get_data().push_back(new_column);
            }
        }

        for (auto& [label, node_matrix] : result)
        {
            node_matrix.preprocess();
        }
        return result;
#elif SEQ_TYPE_AA
        throw std::runtime_error("PhyML for proteins is not supported yet.");
#else
        static_assert(false, """Make sure the sequence type is defined. Supported types:\n"""
                             """SEQ_TYPE_DNA"""
                             """SEQ_TYPE_AA""");
#endif

    }

    raxmlng_reader::raxmlng_reader(const string& file_name) noexcept
        : _file_name{ file_name }
    {}

    proba_matrix raxmlng_reader::read()
    {
        try
        {
            cout << "Loading RAXML-NG results: " + _file_name << "..." << endl;
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            auto matrix = read_matrix();

            cout << "Loaded " << matrix.num_branches() << " matrices of " <<
                 matrix.num_sites() << " rows." << endl;
            cout << "Time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - begin).count() << endl << endl;
            return matrix;
        }
        catch (::io::error::integer_overflow& error)
        {
            throw std::runtime_error("PhyML result parsing error: " + string(error.what()));
        }
    }

    proba_matrix raxmlng_reader::read_matrix()
    {
        proba_matrix result;

#ifdef SEQ_TYPE_DNA
        ::io::CSVReader<5,
            ::io::trim_chars<' '>,
            ::io::no_quote_escape<'\t'>,
            ::io::throw_on_overflow,
            ::io::single_and_empty_line_comment<'.'>> _in(_file_name);

        _in.read_header(::io::ignore_extra_column, "Node", "p_A", "p_C", "p_G", "p_T");

        std::string node_label;
        phylo_kmer::score_type a, c, g, t;
        while (_in.read_row(node_label, a, c, g, t))
        {
            auto new_column = std::vector<phylo_kmer::score_type>{ a, c, g, t };

            /// log-transform the probabilities
            auto log = [](auto value) { return std::log10(value); };
            std::transform(begin(new_column), end(new_column), begin(new_column), log);

            auto& node_matrix = result[node_label];
            node_matrix.set_label(node_label);
            node_matrix.get_data().push_back(new_column);
        }

        for (auto& [label, node_matrix] : result)
        {
            node_matrix.preprocess();
        }

#elif SEQ_TYPE_AA
        ::io::CSVReader<21,
            ::io::trim_chars<' '>,
            ::io::no_quote_escape<'\t'>,
            ::io::throw_on_overflow,
            ::io::single_and_empty_line_comment<'.'>> _in(_file_name);

        _in.read_header(::io::ignore_extra_column, "Node", "p_A", "p_R", "p_N", "p_D", "p_C", "p_Q", "p_E", "p_G",
                        "p_H", "p_I", "p_L", "p_K", "p_M", "p_F", "p_P", "p_S", "p_T", "p_W", "p_Y", "p_V");

        std::string node_label;
        i2l::phylo_kmer::score_type a, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v;
        while (_in.read_row(node_label, a, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v))
        {
            /// the order of acids in the RAxML-ng format is not the same as
            /// in the encoding of RAPPAS and IPK
            auto new_column = std::vector<phylo_kmer::score_type> {
                r, h, k, d, e, s, t, n, q, c, g, p, a, i, l, m, f, w, y, v
            };

            /// log-transform the probabilities
            auto log = [](auto value) { return std::log10(value); };
            std::transform(begin(new_column), end(new_column), begin(new_column), log);

            auto& node_matrix = result[node_label];
            node_matrix.set_label(node_label);
            node_matrix.get_data().push_back(new_column);
        }

        for (auto& [label, node_matrix] : result)
        {
            node_matrix.preprocess();
        }
#else
            static_assert(false, """Make sure the sequence type is defined. Supported types:\n"""
                             """SEQ_TYPE_DNA"""
                             """SEQ_TYPE_AA""");
#endif

        return result;
    }

    /// Figures out which AR software is used by running BINARY_FILE --help
    class ar_guesser
    {
    public:
        ar_guesser(std::string working_directory, std::string binary_file)
            : _working_dir{ std::move(working_directory) }
            , _binary_file{ std::move(binary_file) }
        {
            fs::create_directories(_working_dir);
            _ar_output_file = fs::path{ _working_dir } / "ar_help.log";
        }

        ar::software run()
        {
            try
            {
                bp::system(_binary_file, "--help", bp::std_out > _ar_output_file.string());
                return parse_ar_output();
            }
            catch (const std::system_error& err)
            {
                throw std::runtime_error("Error: Could not run ancestral reconstruction software: " + _binary_file);
            }

        }

    private:
        ar::software parse_ar_output()
        {
            std::ifstream in(_ar_output_file.string());

            if (in.is_open())
            {
                std::string line;
                while (getline(in, line))
                {
                    boost::algorithm::to_lower(line);

                    if (line.find("phyml") != string::npos)
                    {
                        return ar::software::PHYML;
                    }
                    else if (line.find("raxml-ng") != string::npos)
                    {
                        return ar::software::RAXML_NG;
                    }
                }
            }

            throw std::runtime_error("Error: Unsupported ancestral reconstruction software: " + _binary_file);
        }

        std::string _working_dir;
        std::string _binary_file;

        fs::path _ar_output_file;
    };

    std::unique_ptr<ar_reader> make_reader(ar::software software, const string& filename)
    {
        if (software == ar::software::PHYML)
        {
            return std::make_unique<phyml_reader>(filename);
        }
        else if (software == ar::software::RAXML_NG)
        {
            return std::make_unique<raxmlng_reader>(filename);
        }
        else
        {
            throw std::runtime_error("Unsupported ancestral reconstruction output format.");
        }
    }

/*
    /// \brief Reads a "extended_tree_node_mapping.tsv" file produced by the old RAPPAS.
    extended_mapping load_extended_mapping(const string& file_name)
    {
        cout << "Loading a node mapping: " + file_name << endl;
        extended_mapping mapping;

        ::io::CSVReader<2, ::io::trim_chars<' '>, ::io::no_quote_escape<'\t'>> in(file_name);
        in.read_header(::io::ignore_extra_column, "original_id", "extended_name");
        std::string extended_name;
        branch_type original_id = i2l::phylo_kmer::na_branch;
        while (in.read_row(original_id, extended_name))
        {
            mapping[extended_name] = original_id;
        }
        cout << "Loaded " << mapping.size() << " mapped ids." << endl << endl;
        return mapping;
    }*/

    std::optional<i2l::phylo_kmer::branch_type> extract_number(const std::string& s)
    {
        /// take the first sequence of digits from the string
        string output = std::regex_replace(s,
            std::regex("[^0-9]*([0-9]+).*"),
            string("$1")
        );

        /// return it if there was any
        if (output.size() > 0)
        {
            const auto casted = static_cast<i2l::phylo_kmer::branch_type>(std::stoul(output));
            return { casted };
        }
        else
        {
            return std::nullopt;
        }
    }

    /// \brief Returns if the input string is a number
    bool is_number(const std::string& s)
    {
        auto it = s.begin();
        while (it != s.end() && std::isdigit(*it))
        {
            ++it;
        }
        return !s.empty() && it == s.end();
    }
    /*

    /// \brief Reads a "ARtree_id_mapping.tsv" file produced by the old RAPPAS.
    artree_label_mapping load_artree_mapping(const std::string& file_name)
    {
        cout << "Loading a node mapping: " + file_name << endl;
        artree_label_mapping mapping;

        ::io::CSVReader<2, ::io::trim_chars<' '>, ::io::no_quote_escape<'\t'>> in(file_name);
        in.read_header(::io::ignore_extra_column, "extended_label", "ARtree_label");
        std::string extended_label, artree_label;
        while (in.read_row(extended_label, artree_label))
        {
            mapping[extended_label] = artree_label;
        }
        cout << "Loaded " << mapping.size() << " mapped ids." << endl << endl;
        return mapping;
    }
*/

    ar::model parse_model(const std::string& model)
    {
        std::map<std::string, ar::model> map = {
            {"JC", ar::model::JC },
            {"K80", ar::model::K80 },
            {"F81", ar::model::F81 },
            {"HKY", ar::model::HKY },
            {"F84", ar::model::F84 },
            {"TN93ef", ar::model::TN93ef },
            {"TN93", ar::model::TN93 },
            {"K81", ar::model::K81 },
            {"K81uf", ar::model::K81uf },
            {"TPM2", ar::model::TPM2 },
            {"TPM2uf", ar::model::TPM2uf },
            {"TPM3", ar::model::TPM3 },
            {"TPM3uf", ar::model::TPM3uf },
            {"TIM1", ar::model::TIM1 },
            {"TIM1uf", ar::model::TIM1uf },
            {"TIM2", ar::model::TIM2 },
            {"TIM2uf", ar::model::TIM2uf },
            {"TIM3", ar::model::TIM3 },
            {"TIM3uf", ar::model::TIM3uf },
            {"TVMef", ar::model::TVMef },
            {"TVM", ar::model::TVM },
            {"SYM", ar::model::SYM },
            {"GTR", ar::model::GTR },

            // protein models
            {"JTT", ar::model::JTT },
            {"LG", ar::model::LG },
            {"WAG", ar::model::WAG },
        };

        if (const auto it = map.find(model); it != map.end())
        {
            return it->second;
        }
        else
        {
            throw std::runtime_error("Unsupported AR model: " + model);
        }
    }

    std::string model_to_string(ar::model model)
    {
        switch(model)
        {
            case ar::model::JC:
                return "JC";
            case ar::model::K80:
                return "K80";
            case ar::model::F81:
                return "F81";
            case ar::model::HKY:
                return "HKY";
            case ar::model::F84:
                return "F84";
            case ar::model::TN93ef:
                return "TN93ef";
            case ar::model::TN93:
                return "TN93";
            case ar::model::K81:
                return "K81";
            case ar::model::K81uf:
                return "K81uf";
            case ar::model::TPM2:
                return "TPM2";
            case ar::model::TPM2uf:
                return "TPM2uf";
            case ar::model::TPM3:
                return "TPM3";
            case ar::model::TPM3uf:
                return "TPM3uf";
            case ar::model::TIM1:
                return "TIM1";
            case ar::model::TIM1uf:
                return "TIM1uf";
            case ar::model::TIM2:
                return "TIM2";
            case ar::model::TIM2uf:
                return "TIM2uf";
            case ar::model::TIM3:
                return "TIM3";
            case ar::model::TIM3uf:
                return "TIM3uf";
            case ar::model::TVMef:
                return "TVMef";
            case ar::model::TVM:
                return "TVM";
            case ar::model::SYM:
                return "SYM";
            case ar::model::GTR:
                return "GTR";

            case ar::model::JTT:
                return "JTT";
            case ar::model::LG:
                return "LG";
            case ar::model::WAG:
                return "WAG";

        }
        throw std::runtime_error("Internal error: wrong model");
    }

    struct ar_result
    {
        std::string matrix_file;
        std::string tree_file;
    };

    /// Iterates over files in the directory, looks for the first file with the given suffix
    optional<fs::path> find_file_by_suffix(const fs::path& directory, const std::string& suffix)
    {
        auto it = fs::directory_iterator(directory);
        for (const auto& entry : boost::make_iterator_range(it, {}))
        {
            if (fs::is_regular_file(entry) && boost::ends_with(entry.path().string(), suffix))
            {
                return entry;
            }
        }
        return nullopt;
    }

    /// Wrapper for ancestral reconstruction software.
    class ar_wrapper
    {
    public:
        virtual ~ar_wrapper() noexcept = default;

        /// Runs the software. Returns the output files on success
        virtual ar_result run() = 0;
    };

    class phyml_wrapper : public ar_wrapper
    {
    public:
        explicit phyml_wrapper(ar::parameters parameters)
            : _params{ std::move(parameters) }
        {
        }

        ~phyml_wrapper() noexcept override = default;

        ar_result run() override
        {
            fs::path matrix_file;
            fs::path tree_file;

            /// --ar-dir not provided
            if (_params.ar_dir.empty())
            {
                _run();

                matrix_file = _params.alignment_file + "_phyml_ancestral_seq.txt";
                tree_file = _params.alignment_file + "_phyml_ancestral_tree.txt";

                check_file(matrix_file);
                check_file(tree_file);
            }
            /// look in the directory provided by --ar-dir
            else
            {
                const auto ar_dir = _params.ar_dir;
                if (fs::is_directory(_params.ar_dir))
                {
                    if (auto found_matrix = find_file_by_suffix(_params.ar_dir, "_phyml_ancestral_seq.txt"); found_matrix)
                    {
                        matrix_file = *found_matrix;
                    }
                    else
                    {
                        throw std::runtime_error("Could not find \"*_phyml_ancestral_seq.txt\" in"
                                                 "the folder provided by --ar-dir: " + _params.ar_dir);
                    }

                    if (auto found_tree = find_file_by_suffix(_params.ar_dir, "_phyml_ancestral_tree.txt"); found_tree)
                    {
                        tree_file = *found_tree;
                    }
                    else
                    {
                        throw std::runtime_error("Could not find \"*_phyml_ancestral_tree.txt\" in"
                                                 "the folder provided by --ar-dir: " + _params.ar_dir);
                    }
                }
                else
                {
                    throw std::runtime_error("Error! No such directory: " + _params.ar_dir);
                }
            }

            std::cout << "Ancestral reconstruction results have been found: " << std::endl
                      << '\t' << matrix_file.string() << std::endl
                      << '\t' << tree_file.string() << std::endl;
            return { matrix_file.string(), tree_file.string() };
        }

    private:

        void _run()
        {
            int result = bp::system(_params.binary_file,
                                    "--ancestral",
                                    "--no_memory_check",
                                    "-i", _params.alignment_file,
                                    "-u", _params.tree_file,
                                    "-m", model_to_string(_params.ar_model),
                                    "-c", std::to_string(_params.categories),
                                    "-b", std::to_string(0),
                                    "-v", std::to_string(0.0),
                                    "-o", "r",
                                    "-a", std::to_string(_params.alpha),
                                    "-f", "e",
                                    "--leave_duplicates"
            );

            if (result != 0)
            {
                throw std::runtime_error("Error during ancestral reconstruction: exit code "
                                         + std::to_string(result));
            }
        }

        void check_file(const fs::path& file)
        {
            if (!fs::exists(file) || fs::is_empty(file))
            {
                throw std::runtime_error("Error during ancestral reconstruction: could not find "
                                         + file.string());
            }
        }

        ar::parameters _params;

    };

    class raxml_wrapper : public ar_wrapper
    {
    public:
        explicit raxml_wrapper(ar::parameters parameters)
            : _params{std::move(parameters)}
        {
        }

        ~raxml_wrapper() noexcept override = default;

        ar_result run() override
        {
            fs::path matrix_file;
            fs::path tree_file;

            /// --ar-dir not provided
            if (_params.ar_dir.empty())
            {
                _run();

                matrix_file = _params.alignment_file + ".raxml.ancestralProbs";
                tree_file = _params.alignment_file + ".raxml.ancestralTree";

                check_file(matrix_file);
                check_file(tree_file);
            }
                /// look in the directory provided by --ar-dir
            else
            {
                const auto ar_dir = _params.ar_dir;
                if (fs::is_directory(_params.ar_dir))
                {
                    if (auto found_matrix = find_file_by_suffix(_params.ar_dir, ".raxml.ancestralProbs"); found_matrix)
                    {
                        matrix_file = *found_matrix;
                    }
                    else
                    {
                        throw std::runtime_error("Could not find \"*.raxml.ancestralTree\" in"
                                                 "the folder provided by --ar-dir: " + _params.ar_dir);
                    }

                    if (auto found_tree = find_file_by_suffix(_params.ar_dir, ".raxml.ancestralTree"); found_tree)
                    {
                        tree_file = *found_tree;
                    }
                    else
                    {
                        throw std::runtime_error("Could not find \"*.raxml.ancestralTree\" in"
                                                 "the folder provided by --ar-dir: " + _params.ar_dir);
                    }
                }
                else
                {
                    throw std::runtime_error("Error! No such directory: " + _params.ar_dir);
                }
            }

            std::cout << "Ancestral reconstruction results have been found: " << std::endl
                      << '\t' << matrix_file.string() << std::endl
                      << '\t' << tree_file.string() << std::endl;
            return {matrix_file.string(), tree_file.string()};
        }

    private:

        void _run()
        {
            auto process = make_process();
            process.wait();
            auto result = process.exit_code();

            if (result != 0)
            {
                throw std::runtime_error("Error during ancestral reconstruction: exit code "
                                         + std::to_string(result));
            }
        }

        std::string get_model_type(const ar::model& model)
        {
            // return "AA";
            return "DNA";
        }

        bp::child make_process()
        {
            std::vector<std::string> args = {
                "--ancestral",
                "--msa", _params.alignment_file,
                "--tree", _params.tree_file,
                "--threads", _params.num_threads,
                "--precision", "9",
                "--seed", "1",
                "--force", "msa",
                "--redo"
            };

            if (_params.ar_parameters.empty())
            {
                // See: https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model
                args.push_back("--data-type");
                args.push_back(get_model_type(_params.ar_model));

                args.push_back("--model");
                const auto ar_model = std::vector<std::string>{
                    model_to_string(_params.ar_model),
                    "+G", std::to_string(_params.categories),
                    "{", std::to_string(_params.alpha), "}",
                    "+IU{0}",
                    "+FC"
                };
                const auto ar_model_str = boost::algorithm::join(ar_model, "");
                args.push_back(ar_model_str);

                args.push_back("--blopt");
                args.push_back("nr_safe");
                args.push_back("--opt-model");
                args.push_back("on");
                args.push_back("--opt-branches");
                args.push_back("on");
            }
            else
            {
                std::vector<std::string> config_args;
                boost::split(config_args, _params.ar_parameters, boost::is_any_of(" "));
                for (const auto& arg : config_args)
                {
                    args.push_back(arg);
                }
            }
            std::cout << "Running: " << _params.binary_file << " " << boost::algorithm::join(args, " ") << std::endl;
            return bp::child(_params.binary_file, bp::args(args));
        }

        void check_file(const fs::path& file)
        {
            if (!fs::exists(file) || fs::is_empty(file))
            {
                throw std::runtime_error("Error during ancestral reconstruction: could not find "
                                         + file.string());
            }
        }

        ar::parameters _params;

    };

    std::unique_ptr<ar_wrapper> make_ar_wrapper(ar::software software, const ar::parameters& parameters)
    {
        if (software == ar::software::PHYML)
        {
            return std::make_unique<phyml_wrapper>(parameters);
        }
        else if (software == ar::software::RAXML_NG)
        {
            return std::make_unique<raxml_wrapper>(parameters);
        }
        else
        {
            throw std::runtime_error("Unsupported ancestral reconstruction output format.");
        }
    }

    std::pair<ar::software, ar::parameters> make_parameters(const cli::parameters& parameters,
                                                            const std::string& ext_tree_file,
                                                            const std::string& ext_alignment_phylip)
    {
        ar::parameters ar_params;
        ar_params.ar_dir = parameters.ar_dir;
        ar_params.binary_file = parameters.ar_binary_file;
        ar_params.ar_parameters = parameters.ar_parameters;
        ar_params.ar_model = parse_model(parameters.ar_model);
        ar_params.alpha = parameters.ar_alpha;
        ar_params.categories = parameters.ar_categories;
        ar_params.num_threads = std::to_string(parameters.num_threads);
        ar_params.tree_file = ext_tree_file;
        ar_params.alignment_file = ext_alignment_phylip;

        /// Find out what AR software is used
        ar_guesser arguesser(parameters.working_directory, parameters.ar_binary_file);
        const auto ar_software = arguesser.run();

        return { ar_software, ar_params };
    }

    std::tuple<proba_matrix, i2l::phylo_tree> ancestral_reconstruction(ar::software software, const ar::parameters& parameters)
    {
        /// Run ancestral reconstruction
        auto wrapper = make_ar_wrapper(software, parameters);
        const auto& result = wrapper->run();

        /// Parse the probability matrix
        auto reader = make_reader(software, result.matrix_file);
        auto matrix = reader->read();

        /// Read the tree generated by AR software
        auto ar_tree = i2l::io::load_newick(result.tree_file);

        return { std::move(matrix), std::move(ar_tree) };
    }

    size_t count_leaves(const i2l::phylo_tree& tree)
    {
        size_t count = 0;
        for (const auto& node : tree)
        {
            if (node.is_leaf())
            {
                ++count;
            }
        }
        return count;
    }

    ar::mapping map_nodes(const i2l::phylo_tree& extended_tree, const i2l::phylo_tree& ar_tree)
    {
        if (extended_tree.get_node_count() != ar_tree.get_node_count())
        {
            throw std::runtime_error("Error during database construction: extended tree and "
                                     "AR differ in the number of nodes: "
                                     + std::to_string(extended_tree.get_node_count())
                                     + " vs. " + std::to_string(ar_tree.get_node_count()));
        }

        ar::mapping ext_to_ar;

        /// Traverse the AR tree
        using const_iterator = i2l::postorder_tree_iterator<true>;
        for (const auto& ext_node : i2l::visit_subtree<const_iterator>(extended_tree.get_root()))
        {
            if (ext_node.is_root())
            {
                continue;
            }

            if (ext_node.is_leaf())
            {
                /// For leaves the node labels must be the same
                const auto label = ext_node.get_label();
                const auto ar_node_ptr = ar_tree.get_by_label(label);

                if (!ar_node_ptr || *ar_node_ptr == nullptr)
                {
                    throw std::runtime_error("Internal error: could not find a node in the AR tree: label = '"
                                             + label + "'");
                }

                /// Map the the Extended tree leaf to the AR tree leaf
                /// They are actually the same.
                ext_to_ar[label] = (*ar_node_ptr)->get_label();

                /// Map their parents if they exist
                if (!ext_node.is_root())
                {
                    const auto ext_parent_ptr = ext_node.get_parent();
                    const auto ar_parent_ptr = (*ar_node_ptr)->get_parent();

                    ext_to_ar[ext_parent_ptr->get_label()] = ar_parent_ptr->get_label();
                }
            }
            /// Internal node
            else
            {
                const auto ar_node_label = ext_to_ar[ext_node.get_label()];
                const auto ar_node_ptr = ar_tree.get_by_label(ar_node_label);
                if (!ar_node_ptr || *ar_node_ptr == nullptr)
                {
                    throw std::runtime_error("Internal error: could not a find node in the AR tree: label = '"
                                             + ar_node_label + "'");
                }

                /// Map their parent if they exist
                const auto ext_parent_ptr = ext_node.get_parent();
                const auto ar_parent_ptr = (*ar_node_ptr)->get_parent();
                if (ar_parent_ptr && ext_parent_ptr)
                {
                    ext_to_ar[ext_parent_ptr->get_label()] = ar_parent_ptr->get_label();
                }
            }
        }

        /*std::cout << std::endl << "EXT -> AR MAPPING" << std::endl;
        for (const auto& [ext_label, ar_label] : ext_to_ar)
        {
            std::cout << ext_label << " -> " << ar_label << std::endl;
        }*/
        return ext_to_ar;
    }
}
