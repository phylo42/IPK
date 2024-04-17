#include <iostream>
#include <i2l/phylo_kmer_db.h>
#include <i2l/phylo_tree.h>
#include <i2l/newick.h>
#include <i2l/serialization.h>

using db = i2l::phylo_kmer_db;


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " DATABASE" << std::endl;
        return 1;
    }

    const auto db = i2l::load(argv[1]);
    const auto tree = i2l::io::parse_newick(db.tree());

    for (const auto& [kmer, entries] : db)
    {
        std::cout << i2l::decode_kmer(kmer, db.kmer_size()) << std::endl;
#if defined(KEEP_POSITIONS)
        for (const auto& [branch, score, position] : entries)
#else
        for (const auto& [branch, score] : entries)
#endif
        {
            const auto node = *tree.get_by_postorder_id(branch);
            std::cout << "\t" << std::pow(10, score) << "\t" << node->get_preorder_id() << std::endl;
        }
    }

}