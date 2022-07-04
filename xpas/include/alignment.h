#ifndef RAPPAS_BUILD_ALIGNMENT_H
#define RAPPAS_BUILD_ALIGNMENT_H

#include <vector>
#include <string>
#include <memory>
#include <xcl/seq.h>
#include <xcl/seq_record.h>

namespace xcl
{
    class phylo_tree;
}

namespace xpas
{
    using xcl::seq_record;
    class alignment;

    alignment extend_alignment(alignment original_alignment, const xcl::phylo_tree& tree);

    class alignment
    {
        friend alignment extend_alignment(alignment original_alignment, const xcl::phylo_tree& tree);
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

    alignment extend_alignment(alignment original_alignment, const xcl::phylo_tree& tree);


    enum class alignment_format
    {
        FASTA,
        PHYLIP
    };

    void save_alignment(const alignment& alignment, const std::string& filename,
                        alignment_format format);
}


#endif //RAPPAS_BUILD_ALIGNMENT_H
