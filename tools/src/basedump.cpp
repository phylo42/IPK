#include <iostream>
#include <xcl/phylo_kmer_db.h>
#include <xcl/serialization.h>

/// For a given databases, prints k-mers and stored positions
class basedump
{
public:
    basedump(const std::string& filename)
        : db(xcl::load(filename))
    {
    }

    void dump()
    {
        std::cout << "kmer\tbranch\tscore\tpos" << std::endl;
        for (const auto& [kmer, entries] : db)
        {
            for (const auto& [branch, score, pos] : entries)
            {
                std::cout << xcl::decode_kmer(kmer, db.kmer_size()) << "\t" << branch << "\t" << score << "\t" << pos << std::endl;
            }
        }
    }

private:
    xcl::phylo_kmer_db db;
};


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cout << "Usage:\n\t" << argv[0] << " DB_FILE" << std::endl;
        return -1;
    }

    basedump bd(argv[1]);
    bd.dump();
    return 0;
}