#include <iostream>
#include <boost/filesystem.hpp>
#include <xpas/phylo_kmer_db.h>
#include <xpas/serialization.h>

uintmax_t get_num_entries(const xpas::phylo_kmer_db& db)
{
    uintmax_t num_entries = 0;

    for (const auto& [key, entries] : db)
    {
        num_entries += entries.size();
    }
    return num_entries;
}

std::string humanize_size(boost::uintmax_t value)
{
    const auto factors = std::vector<std::string>{ "b", "Kb", "Mb", "Gb", "Tb" };
    int i = 0;
    while (value > 1024)
    {
        value /= 1024;
        i++;
    }
    return std::to_string(value) + " " + factors[i];
}


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " DATABASE_FILENAME" << std::endl;
        return 1;
    }

    const auto filename = std::string{ argv[1] };
    std::cout << "Loading " << filename << "..." << std::endl;
    const auto db = xpas::load(filename);
    const auto num_keys = db.size();
    const auto num_entries = get_num_entries(db);
    const auto filesize = boost::filesystem::file_size(filename);
    const auto estimated_size = num_keys * sizeof(xpas::phylo_kmer::key_type) +
                                sizeof(xpas::pkdb_value) * num_entries;
    std::cout << "Database parameters:" << std::endl
              << "\tSequence type: " << db.sequence_type() << std::endl
              << "\tk: " << db.kmer_size() << std::endl
              << "\tomega: " << db.omega() << std::endl
              << "\tKeep positions: " << (db.positions_loaded() ? "true" : "false") << std::endl
              << "\t# k-mers: " << num_keys << std::endl
              << "\t# k-mer entries: " << num_entries << std::endl
              << "\tDisk size (compressed): " << humanize_size(filesize) << std::endl
              << "\tEstimated RAM size: ~" << humanize_size(estimated_size) << std::endl
              << std::endl;
}