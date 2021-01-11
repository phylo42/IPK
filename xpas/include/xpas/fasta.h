#ifndef XPAS_FASTA_H
#define XPAS_FASTA_H


#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <string>
#include <vector>
#include <xpas/seq_record.h>

namespace xpas::io
{
    namespace impl
    {
        namespace bio = boost::iostreams;

        /// Fasta file iterator, implements batch sequence reading.
        /// WARNING: Returns references to sequences stored in the iterator.
        /// Do not keep these references after usage.
        class fasta_iterator
        {
        public:
            using value_type = const xpas::seq_record&;

            fasta_iterator(const std::string& filename, size_t batch_size, bool clean_sequences=true);
            fasta_iterator(const fasta_iterator&) = delete;
            fasta_iterator(fasta_iterator&&) = default;
            fasta_iterator& operator=(const fasta_iterator&) = delete;
            fasta_iterator& operator=(fasta_iterator&&) = delete;
            ~fasta_iterator() noexcept = default;

            fasta_iterator& operator++();

            bool operator==(const fasta_iterator& rhs) const noexcept;
            bool operator!=(const fasta_iterator& rhs) const noexcept;

            value_type operator*() const noexcept;
        private:
            void _read_batch();

            bio::mapped_file_source _mmap;
            bio::stream<bio::mapped_file_source> _is;

            /// current batch of sequences
            std::vector<seq_record> _seqs;

            size_t _batch_size;
            /// the index of the current sequence in the batch
            size_t _seq_id;
            /// the global index of the current sequence
            size_t _global_seq_id;
            bool _last_batch;

            std::string _header;

            /// Indicates if we need to clean sequences on-the-fly with xpas::clean_sequence
            bool _clean_sequences;
        };
    }

    /// Reads fasta file in batches. Cleans the sequences with clean_sequence
    /// Usage:
    ///     for (const auto& seq: xpas::io::read_fasta(filename)) { ... }
    class read_fasta
    {
        using const_iterator = impl::fasta_iterator;
    public:
        explicit read_fasta(std::string filename, bool clean_sequences=true, size_t batch_size=1024);

        [[nodiscard]]
        const_iterator begin() const;
        [[nodiscard]]
        const_iterator end() const;

    private:
        std::string _filename;
        size_t _batch_size;
        bool _clean_sequences;
    };

    /// Clean an input sequence from gaps. Characters '*', '!', '.' are also skipped
    /// and ignored
    std::string clean_sequence(std::string sequence);
}


/// \brief Outputs a collection of fasta records
std::ostream& operator<<(std::ostream& out, const std::vector<xpas::seq_record>& sequences);

#endif //RAPPAS_CORE_FASTA_H
