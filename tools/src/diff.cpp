#include <iostream>
#include <xpas/phylo_kmer_db.h>
#include <xpas/serialization.h>

using db = xpas::phylo_kmer_db;

class diff
{
public:
    diff(const std::string& filename1, const std::string& filename2)
        : a_filename(filename1), b_filename(filename2),
        a(xpas::load(filename1)), b(xpas::load(filename2))
    {
    }

    template <class T>
    struct check_result
    {
        bool match;
        T a_value;
        T b_value;
    };

    int check()
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

        {
            const auto& [match, va, vb] = check_kmer_size();
            std::cout << "k-mer size:\t" << bool_to_OK(match) << "\t"
                      << va << "\t" << vb << std::endl;
        }

        {
            const auto& [match, va, vb] = check_omega();
            std::cout << "Omega:\t" << bool_to_OK(match) << "\t"
                      << va << "\t" << vb << std::endl;
        }

        {
            const auto& [match, va, vb] = check_tree();
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




private:
    std::string a_filename;
    std::string b_filename;
    db a;
    db b;
};




int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::cout << "Usage:\n\t" << argv[0] << " DB_FILE1 DB_FILE2" << std::endl;
        return -1;
    }

    diff checker(argv[1], argv[2]);
    return checker.check();
}