#include <iostream>
#include <i2l/phylo_kmer_db.h>
#include <i2l/serialization.h>

using db = i2l::phylo_kmer_db;

class diff
{
public:
    diff(const std::string& filename1, const std::string& filename2)
        : a_filename(filename1), b_filename(filename2),
        a(i2l::load(filename1)), b(i2l::load(filename2))
    {
    }

    template <class T>
    struct check_result
    {
        bool match;
        T a_value;
        T b_value;
    };

    int check(bool verbose)
    {
        {
            const auto& [match, va, vb] = check_sequence_type();
            std::cout << "Sequence type:\t" << bool_to_OK(match)
                      << "\t" << va << "\t" << vb << std::endl;
        }

        {
            const auto& [match, va, vb] = check_positions();
            std::cout << "Position support:\t" << bool_to_OK(match) << "\t"
                      << bool_to_str(va) << "\t" << bool_to_str(vb) << std::endl;
        }

        {
            const auto& [match, va, vb] = check_version();
            std::cout << "Protocol version:\t" << bool_to_OK(match) << "\t"
                      << va << "\t" << vb << std::endl;
        }


        const auto& [k_match, a_k, b_k] = check_kmer_size();
        std::cout << "k-mer size:\t" << bool_to_OK(k_match) << "\t"
                  << a_k << "\t" << b_k << std::endl;


        {
            const auto& [match, va, vb] = check_omega();
            std::cout << "Omega:\t" << bool_to_OK(match) << "\t"
                      << va << "\t" << vb << std::endl;

            auto eps = [](float omega, size_t k) { return std::log10(i2l::score_threshold(omega, k)); };
            std::cout << "Threshold:\t" << bool_to_OK(match) << "\t"
                      << eps(va, a_k) << "\t" << eps(vb, b_k) << std::endl;
        }

        {
            const auto& [match, va, vb] = check_tree();
            (void)va; (void) vb;
            std::cout << "Reference tree:\t" << bool_to_OK(match) << "\t"
                      << " " << "\t" << " " << std::endl;
        }

        {
            // TODO: implement it
            std::cout << "Tree index:\t???" << std::endl;
        }

        {
            const auto& [match, va, vb] = check_size();
            std::cout << "Number of k-mers:\t" << bool_to_OK(match) << "\t"
                      << va << "\t" << vb << std::endl;
        }

        {
            const auto& [match, va, vb] = check_num_phylokmers();
            std::cout << "Number of phylo-k-mers:\t" << bool_to_OK(match) << "\t"
                      << va << "\t" << vb << std::endl;
        }

        {
            const auto& [match, diffs] = check_phylo_kmers();
            std::cout << "Phylo-k-mer scores:\t" << bool_to_OK(match) << std::endl;

            if (verbose)
            {
                auto val_or_nan_to_string = [](auto score) { return std::isnan(score) ? "-" : std::to_string(score); };
                (void)val_or_nan_to_string;
                std::cout << "\t\tcode\tk-mer\tbranch\tA score\tB score\n";
                for (const auto& [kmer, branch, a_score, b_score] : diffs)
                {
                    std::cout << "\t\t" << kmer << "\t" << i2l::decode_kmer(kmer, a.kmer_size()) << "\t" << branch << "\t"
                              //<< val_or_nan_to_string(a_score) << "\t"
                              //<< val_or_nan_to_string(b_score) << "\t" << std::endl;
                              << std::pow(10, a_score) << "\t"
                              << std::pow(10, b_score) << "\t" << std::endl;
                }
            }
        }
        return 0;
    }

    static std::string bool_to_OK(bool x)
    {
        return x ? "OK" : "DIFF";
    };

    static std::string bool_to_str(bool x)
    {
        return x ? "true" : "false";
    };

    check_result<std::string_view> check_sequence_type()
    {
        bool match = a.sequence_type() == b.sequence_type();
        return { match, a.sequence_type(), b.sequence_type() };
    }

    check_result<bool> check_positions()
    {
        bool match = a.positions_loaded() == b.positions_loaded();
        return { match, a.positions_loaded(), b.positions_loaded() };
    }

    check_result<unsigned int> check_version()
    {
        bool match = a.version() == b.version();
        return { match, a.version(), b.version() };
    }

    check_result<size_t> check_kmer_size()
    {
        bool match = a.kmer_size() == b.kmer_size();
        return { match, a.kmer_size(), b.kmer_size() };
    }

    check_result<double> check_omega()
    {
        bool match = a.omega() == b.omega();
        return { match, a.omega(), b.omega() };
    }

    check_result<std::string_view> check_tree()
    {
        bool match = a.tree() == b.tree();
        return { match, a.tree(), b.tree() };
    }

    check_result<size_t> check_size()
    {
        bool match = a.size() == b.size();
        return { match, a.size(), b.size() };
    }

    static size_t get_num_phylokmers(const db& x)
    {
        size_t i = 0;
        for (const auto& [kmer, entries] : x)
        {
            (void)kmer;
            i += entries.size();
        }
        return i;
    }

    check_result<size_t> check_num_phylokmers()
    {
        size_t a_size = get_num_phylokmers(a);
        size_t b_size = get_num_phylokmers(b);
        return { a_size == b_size, a_size, b_size };
    }

    using pk_map = std::unordered_map<i2l::phylo_kmer::branch_type, i2l::phylo_kmer::score_type>;

    template<class T>
    pk_map to_map(const T& x)
    {
        pk_map va;
        for (const auto& [a_branch, a_score] : x)
        {
            va[a_branch] = a_score;
        }
        return va;
    }

    struct pk_diff
    {
        i2l::phylo_kmer::key_type kmer;
        i2l::phylo_kmer::branch_type branch;
        i2l::phylo_kmer::score_type a_value;
        i2l::phylo_kmer::score_type b_value;
    };

    std::tuple<bool, std::vector<pk_diff>> check_phylo_kmers()
    {
        const double EPS = 1e-6;

        bool match = true;
        std::vector<pk_diff> diffs;

        /// Check the k-mers of A
        for (const auto& [kmer, a_entries] : a)
        {
            const auto va = to_map(a_entries);

            /// Search for the k-mer in B
            if (auto b_entries = b.search(kmer); b_entries)
            {
                const auto vb = to_map(*b_entries);

                for (const auto& [a_branch, a_score] : a_entries)
                {
                    /// (kmer, branch) is scored in both A and B, check the values
                    if (const auto& it = vb.find(a_branch); it != va.end())
                    {
                        const auto b_score = it->second;
                        bool score_match = std::fabs(a_score - b_score) < EPS;

                        /// scored differently
                        if (!score_match)
                        {
                            diffs.push_back({kmer, a_branch, a_score, b_score});
                        }

                        match &= score_match;
                    }
                    /// the k-mer exists in both A and B, but in A there is a branch
                    /// that is not scored in B
                    else
                    {
                        diffs.push_back({kmer, a_branch, a_score, i2l::phylo_kmer::na_score});
                        match = false;
                    }
                }

                /// Now, check the branches and scores for this k-mer in B
                for (const auto& [b_branch, b_score] : *b_entries)
                {
                    /// the k-mer exists in both A and B, but in B there is a branch
                    /// that is not scored in A
                    if (const auto& it = va.find(b_branch); it == va.end())
                    {
                        diffs.push_back({ kmer, b_branch, i2l::phylo_kmer::na_score, b_score });
                        match = false;
                    }
                }

            }
            /// this k-mer does not exist in B, report all branches
            else
            {
                match = false;

                for (const auto& [a_branch, a_score] : a_entries)
                {
                    diffs.push_back({ kmer, a_branch, a_score, i2l::phylo_kmer::na_score });
                }
            }
        }

        /// We also need to find k-mers in B that do not exist in A.
        /// Those are not covered by the loop above
        for (const auto& [kmer, b_entries] : b)
        {
            auto a_entries = a.search(kmer);
            if (!a_entries)
            {
                /// this k-mer does not exist in A, report all branches
                match = false;

                for (const auto& [b_branch, b_score] : b_entries)
                {
                    diffs.push_back({ kmer, b_branch, i2l::phylo_kmer::na_score, b_score });
                }

            }
        }
        return { match, diffs };
    }

private:
    std::string a_filename;
    std::string b_filename;
    db a;
    db b;
};


int main(int argc, char** argv)
{
    if (argc != 4)
    {
        std::cout << "Usage:\n\t" << argv[0] << "[0/1 VERBOSE] DB_FILE1 DB_FILE2" << std::endl;
        return -1;
    }
    bool verbose = std::string(argv[1]) == "1";

    diff checker(argv[2], argv[3]);
    return checker.check(verbose);
}