#ifndef IPK_PROBA_MATRIX_H
#define IPK_PROBA_MATRIX_H

#include <unordered_map>
#include <memory>
#include <i2l/phylo_kmer.h>
#include "ar.h"
#include "window.h"

namespace ipk
{
    /// \brief A posterior probabilities matrix class.
    /// \details A matrix class for storing posterior probabilities, given by the ancestral reconstruction
    /// algorithm. So, this is a matrix of size [#branch_nodes x #sites x #variants], where:
    /// - #branch_nodes is the number of non-leaf nodes of input tree
    /// - #sites is the size of input alignment,
    /// - #variants is the alphabet size.
    class proba_matrix final
    {
    public:
        using branch_type = i2l::phylo_kmer::branch_type;
        static const branch_type NOT_A_LABEL = std::numeric_limits<branch_type>::max();

        /// to map node labels into corresponding matrices
        using storage = std::unordered_map<std::string, matrix>;
        using iterator = typename storage::iterator;
        using const_iterator = typename storage::const_iterator;
        using mapped_type = storage::mapped_type;

        proba_matrix(std::unique_ptr<ipk::ar::reader> reader);
        proba_matrix(const proba_matrix&) = delete;
        proba_matrix(proba_matrix&& other) = default;
        proba_matrix& operator=(const proba_matrix&) = delete;
        proba_matrix& operator=(proba_matrix&&) = delete;
        ~proba_matrix() = default;

        /// capacity
        [[nodiscard]]
        size_t num_branches() const;

        [[nodiscard]]
        size_t num_sites() const;

        // Lookup
        [[nodiscard]]
        mapped_type& operator[](const std::string& ar_label);

        [[nodiscard]]
        const mapped_type& at(const std::string& ar_label) const;

        [[nodiscard]]
        iterator find(const std::string& ar_label);

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
        storage _data;
        std::unique_ptr<ipk::ar::reader> _reader;
    };
}

#endif