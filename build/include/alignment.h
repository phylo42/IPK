#ifndef RAPPAS_BUILD_ALIGNMENT_H
#define RAPPAS_BUILD_ALIGNMENT_H

#include <vector>
#include <string>
#include <memory>
#include <xpas/seq.h>
#include "fasta.h"

namespace rappas
{
    class alignment
    {
    public:
        typedef std::vector<xpas::io::fasta>::iterator iterator;
        typedef std::vector<xpas::io::fasta>::const_iterator const_iterator;
    public:
        /// ctors
        /// The vector is passed by value and moved
        explicit alignment(std::vector<xpas::io::fasta> sequences);
        alignment(const alignment&) = delete;
        alignment(alignment&&) noexcept = default;
        ~alignment() noexcept = default;
        alignment& operator=(const alignment&) = delete;
        alignment& operator=(alignment&&) noexcept = default;

        /// capacity
        [[nodiscard]]
        size_t width() const noexcept;
        [[nodiscard]]
        size_t height() const noexcept;

        /// iterators
        [[nodiscard]]
        iterator begin();
        [[nodiscard]]
        iterator end();
        [[nodiscard]]
        const_iterator begin() const;
        [[nodiscard]]
        const_iterator end() const;
    private:
        std::vector<xpas::io::fasta> _sequences;

        size_t _width;
        size_t _height;
    };

    /// \brief Read and preprocess the reference alignment.
    /// Performs several check on the sequence content and filters out
    /// columns that have gap ratio >= reduction_ratio
    alignment preprocess_alignment(const std::string& working_dir,
                                   const std::string& alignment_file,
                                   double reduction_ratio,
                                   bool no_reduction);
}



#endif //RAPPAS_BUILD_ALIGNMENT_H
