#include <iostream>
#include <numeric>
#include <cmath>
#include <string>
#include <optional>
#include <regex>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/process.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <csv-parser/csv.h>
#include <xpas/newick.h>
#include <xpas/phylo_tree.h>
#include "ar.h"
#include "row.h"
#include "proba_matrix.h"
#include "command_line.h"

using std::string;
using std::cout, std::endl;
namespace bp = boost::process;
namespace fs = boost::filesystem;

namespace xpas::ar
{
    /// \brief Reads ancestral reconstruction output
    class ar_reader
    {
    public:
        virtual ~ar_reader() noexcept = default;
        virtual xpas::proba_matrix read() = 0;
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

        proba_matrix matrix;
        std::ifstream infile(_file_name);

        throw std::runtime_error("PhyML is not supported in this version of XPAS.");

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
                xpas::phylo_kmer::score_type a, c, g, t;

                std::istringstream iss(line);
                if (!(iss >> site >> node_label >> a >> c >> g >> t))
                {
                    throw std::runtime_error("Parsing error: could not parse the line " + line);
                }

                /// log-transform the probabilities
                auto new_row = row_type { { { a, 0 }, { c, 1 }, { g, 2 }, { t, 3 } } };
                auto log = [](const proba_pair& p) { return proba_pair{ std::log10(p.score), p.index }; };
                std::transform(begin(new_row), end(new_row), begin(new_row), log);

                // sort them
                auto compare = [](const proba_pair& p1, const proba_pair& p2) { return p1.score > p2.score; };
                std::sort(begin(new_row), end(new_row), compare);

                /// insert
                auto it = matrix.find(node_label);
                if (it != std::end(matrix))
                {
                    //it->second.push_back(new_row);
                }
                else
                {
                    /// WARNING: it is cheap to copy new_row here in the case of DNA.
                    /// It is probably makes sense to move it in the case of amino acids though.
                    //matrix[node_label] = proba_matrix::mapped_type{ node_label, { new_row } };
                }
            }
        }
        return matrix;
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
        proba_matrix matrix;

#ifdef SEQ_TYPE_DNA
        ::io::CSVReader<5,
            ::io::trim_chars<' '>,
            ::io::no_quote_escape<'\t'>,
            ::io::throw_on_overflow,
            ::io::single_and_empty_line_comment<'.'>> _in(_file_name);

        _in.read_header(::io::ignore_extra_column, "Node", "p_A", "p_C", "p_G", "p_T");

        proba_matrix result;

        std::string node_label;
        phylo_kmer::score_type a, c, g, t;
        while (_in.read_row(node_label, a, c, g, t))
        {
            auto new_column = std::vector<phylo_kmer::score_type>{ a, c, g, t };

            /// log-transform the probabilities
            auto log = [](auto value) { return std::log10(value); };
            std::transform(begin(new_column), end(new_column), begin(new_column), log);

            result[node_label].get_data().push_back(new_column);
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
        xpas::phylo_kmer::score_type a, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v;
        while (_in.read_row(node_label, a, r, n, d, c, q, e, g, h, i, l, k, m, f, p, s, t, w, y, v))
        {
            /// log-transform the probabilities
            auto new_row = row_type {
                { { a, *xpas::seq_traits::key_to_code('a') },
                    { r, *xpas::seq_traits::key_to_code('r') },
                    { n, *xpas::seq_traits::key_to_code('n')},
                    { d, *xpas::seq_traits::key_to_code('d') },
                    { c, *xpas::seq_traits::key_to_code('c') },
                    { q, *xpas::seq_traits::key_to_code('q') },
                    { e, *xpas::seq_traits::key_to_code('e') },
                    { g, *xpas::seq_traits::key_to_code('g') },
                    { h, *xpas::seq_traits::key_to_code('h') },
                    { i, *xpas::seq_traits::key_to_code('i') },
                    { l, *xpas::seq_traits::key_to_code('l') },
                    { k, *xpas::seq_traits::key_to_code('k') },
                    { m, *xpas::seq_traits::key_to_code('m') },
                    { f, *xpas::seq_traits::key_to_code('f') },
                    { p, *xpas::seq_traits::key_to_code('p') },
                    { s, *xpas::seq_traits::key_to_code('s') },
                    { t, *xpas::seq_traits::key_to_code('t') },
                    { w, *xpas::seq_traits::key_to_code('w') },
                    { y, *xpas::seq_traits::key_to_code('y') },
                    { v, *xpas::seq_traits::key_to_code('v') } }
            };
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
        branch_type original_id = xpas::phylo_kmer::na_branch;
        while (in.read_row(original_id, extended_name))
        {
            mapping[extended_name] = original_id;
        }
        cout << "Loaded " << mapping.size() << " mapped ids." << endl << endl;
        return mapping;
    }*/

    std::optional<xpas::phylo_kmer::branch_type> extract_number(const std::string& s)
    {
        /// take the first sequence of digits from the string
        string output = std::regex_replace(s,
            std::regex("[^0-9]*([0-9]+).*"),
            string("$1")
        );

        /// return it if there was any
        if (output.size() > 0)
        {
            const auto casted = static_cast<xpas::phylo_kmer::branch_type>(std::stoul(output));
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
            {"JC69", ar::model::JC69 },
            {"K80", ar::model::K80 },
            {"F81", ar::model::F81 },
            {"F84", ar::model::F84 },
            {"HKY85", ar::model::HKY85 },
            {"TN93",ar::model::TN93 },
            {"GTR", ar::model::GTR },
            {"LG", ar::model::LG },
            {"WAG", ar::model::WAG },
            {"JTT", ar::model::JTT },
            {"DAYHOFF", ar::model::DAYHOFF },
            {"DCMUT", ar::model::DCMUT },
            {"CPREV", ar::model::CPREV },
            {"MTMAM", ar::model::MTMAM },
            {"MTREV", ar::model::MTREV },
            {"MTART", ar::model::MTART }
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
            case ar::model::JC69:
                return "JC69";
            case ar::model::K80:
                return "K80";
            case ar::model::F81:
                return "F81";
            case ar::model::F84:
                return "F84";
            case ar::model::HKY85:
                return "HKY85";
            case ar::model::TN93:
                return "TN93";
            case ar::model::GTR:
                return "GTR";
            case ar::model::LG:
                return "LG";
            case ar::model::WAG:
                return "WAG";
            case ar::model::JTT:
                return "JTT";
            case ar::model::DAYHOFF:
                return "DAYHOFF";
            case ar::model::DCMUT:
                return "DCMUT";
            case ar::model::CPREV:
                return "CPREV";
            case ar::model::MTMAM:
                return "MTMAM";
            case ar::model::MTREV:
                return "MTREV";
            case ar::model::MTART:
                return "MTART";
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
            : _params{ std::move(parameters) }
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
            return { matrix_file.string(), tree_file.string() };
        }

    private:

        void _run()
        {
            throw std::runtime_error("RAXML-NG is not supported yet.");
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

    std::tuple<proba_matrix, phylo_tree> ancestral_reconstruction(ar::software software, const ar::parameters& parameters)
    {
        /// Run ancestral reconstruction
        auto wrapper = make_ar_wrapper(software, parameters);
        const auto& result = wrapper->run();

        /// Parse the probability matrix
        auto reader = make_reader(software, result.matrix_file);
        auto matrix = reader->read();

        /// Read the tree generated by AR software
        auto ar_tree = xpas::io::load_newick(result.tree_file);

        return { std::move(matrix), std::move(ar_tree) };
    }

    size_t count_leaves(const phylo_tree& tree)
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

    ar::mapping map_nodes(const phylo_tree& extended_tree, const phylo_tree& ar_tree)
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
        using const_iterator = postorder_tree_iterator<true>;
        for (const auto& ext_node : visit_subtree<const_iterator>(extended_tree.get_root()))
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