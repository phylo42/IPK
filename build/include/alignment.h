#ifndef RAPPAS_BUILD_ALIGNMENT_H
#define RAPPAS_BUILD_ALIGNMENT_H

#include <vector>
#include <string>
#include <memory>
#include <xpas/seq.h>

namespace xpas
{
    class alignment;
    class phylo_tree;
    alignment extend_alignment(alignment original_alignment, const phylo_tree& tree);

    class seq_record
    {
    public:
        seq_record() = delete;
        seq_record(std::string header, std::string sequence) noexcept;
        seq_record(seq_record&& other) noexcept = default;
        seq_record(const seq_record&) = default;
        ~seq_record() = default;

        [[nodiscard]]
        std::string_view header() const noexcept;
        [[nodiscard]]
        std::string_view sequence() const noexcept;

        bool operator==(const seq_record& rhs) const noexcept;
    private:
        std::string _header;
        std::string _sequence;
    };

    class alignment
    {
        friend alignment extend_alignment(alignment original_alignment, const phylo_tree& tree);
    public:
        using iterator = std::vector<seq_record>::iterator;
        using const_iterator = std::vector<seq_record>::const_iterator;

    public:
        /// ctors
        /// The vector is passed by value and moved
        explicit alignment(std::vector<seq_record> sequences);
        alignment(const alignment&) = default;
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
        std::vector<seq_record> _sequences;

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

    class phylo_tree;

    alignment extend_alignment(alignment original_alignment, const phylo_tree& tree);


    enum class alignment_format
    {
        FASTA,
        PHYLIP
    };

    void save_alignment(const alignment& alignment, const std::string& filename,
                        alignment_format format);
}


#endif //RAPPAS_BUILD_ALIGNMENT_H
