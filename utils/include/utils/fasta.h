#ifndef RAPPAS_CORE_FASTA_H
#define RAPPAS_CORE_FASTA_H

#include <string>
#include <vector>

namespace utils
{
    class fasta
    {
    public:
        fasta() = delete;
        fasta(std::string&& header, std::string&& sequence) noexcept;
        fasta(fasta&& other) noexcept = default;
        fasta(const fasta&) = delete;
        ~fasta() = default;

        std::string get_header() const noexcept;
        std::string get_sequence() const noexcept;
    private:
        std::string _header;
        std::string _sequence;
    };

    std::vector<fasta> read_fasta(const std::string& filename);

}

#endif //RAPPAS_CORE_FASTA_H
