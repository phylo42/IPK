#include <iostream>
#include <fstream>
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

namespace bp = boost::process;
namespace fs = boost::filesystem;

namespace ipk::ar
{
    /// \brief Reads a PhyML output into a matrix.
    class phyml_reader : public reader
    {
    public:
        phyml_reader(const std::string& file_name) noexcept;
        phyml_reader(const phyml_reader&) = delete;
        phyml_reader(phyml_reader&&) = delete;
        phyml_reader& operator=(const phyml_reader&) = delete;
        phyml_reader& operator=(phyml_reader&&) = delete;
        ~phyml_reader() noexcept override = default;

        ipk::matrix read_node(const std::string& node_label) override;

    private:
        proba_matrix read_matrix();

        std::string _file_name;
    };

    /// \brief Reads RAXML-NG output into a matrix.
    class raxmlng_reader : public reader
    {
    public:
        raxmlng_reader(std::string file_name) noexcept;
        raxmlng_reader(const phyml_reader&) = delete;
        raxmlng_reader(raxmlng_reader&&) = delete;
        raxmlng_reader& operator=(const raxmlng_reader&) = delete;
        raxmlng_reader& operator=(raxmlng_reader&&) = delete;
        ~raxmlng_reader() noexcept override = default;

        ipk::matrix read_node(const std::string& node_label) override;

    private:
        void build_index();
        proba_matrix read_matrix();

        std::string _file_name;
        std::ifstream _file_stream;

        /// Index for Node -> position where the matrix for this node starts
        std::unordered_map<std::string, std::streampos> _index;
    };

    phyml_reader::phyml_reader(const std::string& file_name) noexcept
        : _file_name{ file_name }
    {}


    ipk::matrix phyml_reader::read_node(const std::string& node_label)
    {
        (void)node_label;
        throw std::runtime_error("PhyML is not supported in this version.");
    }
/*
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

                auto new_column = std::array<phylo_kmer::score_type, 4>{ a, c, g, t };

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
            (void)label;
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

    }*/

    raxmlng_reader::raxmlng_reader(std::string file_name) noexcept
        : _file_name{ std::move(file_name) }
    {
        build_index();
    }

    void raxmlng_reader::build_index()
    {
        std::cout << "Indexing " <<  _file_name << "..." << std::endl;
        _file_stream.open(_file_name);

        std::string line;
        std::string current_node;

        /// Skip the header
        std::getline(_file_stream, line);

        /// Remember the stream position right before reading
        auto last_pos = _file_stream.tellg();

        /// A lambda to get the node label from a line
        auto get_label = [](const std::string& line) {
            auto pos = line.find('\t');
            return line.substr(0, pos);
        };

        while (std::getline(_file_stream, line))
        {
            auto node_label = get_label(line);

            /// If read a new label, index the position before
            /// the current line as the beginning of the node block
            if (node_label != current_node)
            {
                _index[node_label] = last_pos;
                current_node = std::move(node_label);
            }

            last_pos = _file_stream.tellg();
        }

        /// Process the last node
        const auto node_label = get_label(line);
        _index[node_label] = last_pos;
    }

    /// The type for the csv reader for a given sequence type
    template<typename SeqType>
    using csv_reader = typename ::io::CSVReader<
        /// Number of columns: Node + Site + State + every char in the alphabet
        3 + i2l::seq_traits_impl<SeqType>::alphabet_size,
        ::io::trim_chars<' '>,
        ::io::no_quote_escape<'\t'>,
        ::io::throw_on_overflow,
        ::io::single_and_empty_line_comment<'.'>>;

    ipk::matrix raxmlng_reader::read_node(const std::string& current_node)
    {
        ipk::matrix matrix;

        /// Number of columns: Node + Site + State + every char in the alphabet
        constexpr size_t num_columns = 3 + i2l::seq_traits::alphabet_size;

        /// The stream position of the node matrix in the file
        const auto pos = _index[current_node];

        /// Make a csv-reader located at the node matrix position
        std::ifstream file_stream(_file_name);
        file_stream.seekg(pos);
        ::io::CSVReader<num_columns,
            ::io::trim_chars<' '>,
            ::io::no_quote_escape<'\t'>,
            ::io::throw_on_overflow,
            ::io::single_and_empty_line_comment<'.'>> _in(_file_name, file_stream);

        bool started = false;
        std::string node_label, site, state;
#if defined(SEQ_TYPE_DNA)
        phylo_kmer::score_type a, c, g, t;
        while (_in.read_row(node_label, site, state, a, c, g, t))
        {
            auto new_column = std::array<phylo_kmer::score_type, 4>{ a, c, g, t };
#elif defined(SEQ_TYPE_AA)
        i2l::phylo_kmer::score_type a, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v;
        while (_in.read_row(node_label, site, state, a, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v))
        {
            /// the order of acids in the RAxML-ng format is not the same as
            /// in the encoding of RAPPAS and IPK
            auto new_column = std::array<phylo_kmer::score_type, 20>{
                r, h, k, d, e, s, t, n, q, c, g, p, a, i, l, m, f, w, y, v
            };

#else
        static_assert(false, """Make sure the sequence type is defined. Supported types:\n"""
                     """SEQ_TYPE_DNA"""
                     """SEQ_TYPE_AA""");
#endif
            (void)site;
            (void)state;
            if (node_label != current_node)
            {
                /// If we did not read anything, that's an error
                [[unlikely]]
                if (!started)
                {
                    throw std::runtime_error("Error while AR indexing: wrong position for node " + current_node);
                }
                /// Finished reading the matrix
                break;
            }

            started = true;

            /// log-transform the probabilities
            auto log = [](auto value) { return std::log10(value); };
            std::transform(begin(new_column), end(new_column), begin(new_column), log);
            matrix.get_data().push_back(new_column);
        }

        if (!started)
        {
            throw std::runtime_error("Could not read the AR matrix for the node " + current_node);
        }
        matrix.set_label(current_node);
        matrix.preprocess();
        return matrix;
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

                    if (line.find("phyml") != std::string::npos)
                    {
                        return ar::software::PHYML;
                    }
                    else if (line.find("raxml-ng") != std::string::npos)
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

    std::unique_ptr<reader> make_reader(ar::software software, const std::string& filename)
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
                throw std::runtime_error("Error during ancestral reconstruction: exit code"
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
                throw std::runtime_error("Error during ancestral reconstruction: exit code"
                                         + std::to_string(result));
            }
        }

        bp::child make_process()
        {
            std::vector<std::string> args = {
                "--ancestral",
                "--msa", _params.alignment_file,
                "--tree", _params.tree_file,
                "--threads", "1", //_params.threads,
                "--precision", "9",
                "--seed", "1",
                "--force", "msa",
                "--redo"
            };

            if (_params.ar_parameters.empty())
            {
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

        auto reader = make_reader(software, result.matrix_file);
        /// Create the wrapper for the AR results
        auto matrix = proba_matrix(std::move(reader));


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

        /// Traverse the extended and AR tree at the same time
        using const_iterator = i2l::postorder_tree_iterator<true>;
        auto visit_ext = i2l::visit_subtree<const_iterator>(extended_tree.get_root());
        auto ext_it = visit_ext.begin();
        auto visit_ar = i2l::visit_subtree<const_iterator>(ar_tree.get_root());
        auto ar_it = visit_ar.begin();
        while (ext_it != visit_ext.end())
        {
            /// Inner nodes usually have no labels. However, we do not
            /// need them in the mapping as it is only for ghost nodes anyway
            if (ext_it->get_label().empty())
            {
                ++ext_it;
                ++ar_it;
                continue;
            }

            ext_to_ar[ext_it->get_label()] = ar_it->get_label();
            ++ext_it;
            ++ar_it;

            if (ext_it->is_root())
            {
                assert(ar_it->is_root());
            }
        }

        assert (ar_it == visit_ar.end());
        return ext_to_ar;
    }
}