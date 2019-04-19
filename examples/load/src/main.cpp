#include <iostream>
#include <core/phylo_kmer_db.h>
#include <core/serialization.h>

core::phylo_kmer_db create_db()
{
    core::phylo_kmer_db db;

    /// branch 0
    db.put(0, 0, 0.0f);
    db.put(1, 0, 0.1f);
    db.put(2, 0, 0.2f);

    /// branch 1
    db.put(1, 1, 0.1f);
    db.put(2, 1, 0.2f);
    db.put(3, 1, 0.3f);

    /// branch 2
    db.put(2, 1, 0.2f);
    db.put(3, 1, 0.3f);
    db.put(4, 1, 0.4f);

    return db;
}

std::ostream& operator<<(std::ostream& out, const core::phylo_kmer_db& db)
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
    const auto filename = "/tmp/rappas_load_example.dat";
    core::save(filename, create_db());

    const auto db = core::load(filename);
    std::cout << db;
}