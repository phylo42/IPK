#include <iostream>
#include <xpas/phylo_kmer_db.h>
#include <xpas/version.h>

xpas::phylo_kmer_db create_db()
{
    const size_t kmer_size = 3;
    const xpas::phylo_kmer::score_type omega = 1.0;
    const std::string tree;

    xpas::phylo_kmer_db db { kmer_size, omega, xpas::seq_type::name, tree };

    /// branch 0
    db.unsafe_insert(0, { 0, 0.00f });
    db.unsafe_insert(1, { 0, 0.10f });
    db.unsafe_insert(2, { 0, 0.20f });

    /// branch 1
    db.unsafe_insert(1, { 1, 0.11f });
    db.unsafe_insert(2, { 1, 0.21f });
    db.unsafe_insert(3, { 1, 0.31f });

    /// branch 2
    db.unsafe_insert(2, { 2, 0.22f });
    db.unsafe_insert(3, { 2, 0.32f });
    db.unsafe_insert(4, { 2, 0.42f });

    return db;
}

std::ostream& operator<<(std::ostream& out, const xpas::phylo_kmer_db& db)
{
    for (const auto& [key, entries] : db)
    {
        out << key << ":\n";
        for (const auto& [branch, score] : entries)
        {
            out << '\t' << branch << ": " << score << '\n';
        }
    }
    return out;
}

int main()
{
    std::cout << "xpas version: " << xpas::version::as_string() << std::endl;

    const auto db = create_db();
    std::cout << "Total number of keys: " << db.size() << std::endl;
    std::cout << "K-mer size: " << db.kmer_size() << std::endl;
    std::cout << db;
}