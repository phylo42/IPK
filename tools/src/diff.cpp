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

    void check()
    {
        check_sequence_type();
        check_positions();
    }

    void check_sequence_type()
    {
        if (a.sequence_type() != b.sequence_type())
        {
            throw std::runtime_error("Sequence types are different: " + std::string(a.sequence_type())
                                     + ", " + std::string(b.sequence_type()));
        }
    }

    void check_positions()
    {
        auto bool_to_str = [](bool x) { return x ? "true" : "false"; };
        if (a.positions_loaded() != b.positions_loaded())
        {
            std::cerr << "One DB has positions while the other has not.\n" <<
                      "\t" << a_filename << ": " << bool_to_str(a.positions_loaded()) << std::endl <<
                      "\t" << b_filename << ": " << bool_to_str(b.positions_loaded()) << std::endl;
        }
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
        std::cout << "Usage:\n\t" << argv[0] << "DB_FILE1 DB_FILE2" << std::endl;
        return -1;
    }

    diff checker(argv[1], argv[2]);
    checker.check();
}