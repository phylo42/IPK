#include <iostream>
#include <xpas/phylo_kmer_db.h>

xpas::phylo_kmer_db create_db()
{
    using namespace xpas;

    const size_t kmer_size = 3;
    const xpas::phylo_kmer::score_type omega = 1.0;
    const std::string tree;

    xpas::phylo_kmer_db db { kmer_size, omega, tree };

    /// branch 0
    db.unsafe_insert(0, make_pkdb_value<phylo_kmer>(0, 0.00f, 0));
    db.unsafe_insert(1, make_pkdb_value<phylo_kmer>(0, 0.10f, 10));
    db.unsafe_insert(2, make_pkdb_value<phylo_kmer>(0, 0.20f, 20));

    /// branch 1
    db.unsafe_insert(1, make_pkdb_value<phylo_kmer>(1, 0.11f, 11));
    db.unsafe_insert(2, make_pkdb_value<phylo_kmer>(1, 0.21f, 21));
    db.unsafe_insert(3, make_pkdb_value<phylo_kmer>(1, 0.31f, 31));

    /// branch 2
    db.unsafe_insert(2, make_pkdb_value<phylo_kmer>(2, 0.22f, 22));
    db.unsafe_insert(3, make_pkdb_value<phylo_kmer>(2, 0.32f, 32));
    db.unsafe_insert(4, make_pkdb_value<phylo_kmer>(2, 0.42f, 42));

    return db;
}

#ifdef KEEP_POSITIONS
void search(const xpas::phylo_kmer_db& db, xpas::phylo_kmer_db::key_type key)
{
    if (auto entries = db.search(key); entries)
    {
        std::cout << "Found " << key << ":\n";
        for (const auto& [branch, score, position] : *entries)
        {
            std::cout << "\tbranch " << branch << ": " << score << " at " << position << '\n';
        }
    }
    else
    {
        std::cout << "Key " << key << " not found.\n";
    }
}
#else
void search(const xpas::phylo_kmer_db& db, xpas::phylo_kmer_db::key_type key)
{
    if (auto entries = db.search(key); entries)
    {
        std::cout << "Found " << key << ":\n";
        for (const auto& [branch, score] : *entries)
        {
            std::cout << "\tbranch " << branch << ": " << score << '\n';
        }
    }
    else
    {
        std::cout << "Key " << key << " not found.\n";
    }
}
#endif
int main()
{
    const auto db = create_db();
    search(db, 0);
    search(db, 2);
    search(db, 42);
}