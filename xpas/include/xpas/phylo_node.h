#ifndef XPAS_BUILD_PHYLO_NODE_H
#define XPAS_BUILD_PHYLO_NODE_H

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <type_traits>

namespace xpas
{
    class phylo_tree;
    class phylo_node;
    class tree_extender;
}

namespace xpas::io
{
    class newick_parser;
    xpas::phylo_tree parse_newick(std::string_view);
}


namespace xpas
{
    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief A node of a phylogenetic tree.
    class phylo_node
    {
        friend xpas::io::newick_parser;
        friend xpas::phylo_tree;
        friend xpas::phylo_tree xpas::io::parse_newick(std::string_view);
        friend xpas::tree_extender;
    public:
        /// Member types

        /// \brief Node post-/pre-order id type
        using id_type = int;

        /// \brief Branch length type
        using branch_length_type = double;

        phylo_node();
        phylo_node(std::string label, branch_length_type branch_length, phylo_node* parent);
        phylo_node(const phylo_node& other) = delete;
        phylo_node& operator=(const phylo_node&) = delete;
        ~phylo_node() noexcept;

        /// WARNING: this operator only checks for the id and label fields
        bool operator==(const phylo_node& rhs) const noexcept;
        bool operator!=(const phylo_node& rhs) const noexcept;

        [[nodiscard]]
        std::string get_label() const noexcept;

        void set_label(const std::string& label);

        [[nodiscard]]
        phylo_node* get_parent() const noexcept;

        [[nodiscard]]
        id_type get_preorder_id() const noexcept;

        [[nodiscard]]
        id_type get_postorder_id() const noexcept;

        [[nodiscard]]
        branch_length_type get_branch_length() const noexcept;

        void set_branch_length(branch_length_type length);

        [[nodiscard]]
        const std::vector<phylo_node*>& get_children() const;

        /// Clean node and fill with the default values. Used in the default constructor
        void clean();

        void add_child(phylo_node* node);

        void remove_child(phylo_node* node);

    private:
        id_type _preorder_id;
        id_type _postorder_id;

        std::string _label;
        branch_length_type _branch_length;

        std::vector<phylo_node*> _children;
        phylo_node* _parent;
    };

}
#endif //XPAS_BUILD_PHYLO_NODE_H
