#include "extended_tree.h"
#include <xpas/newick.h>

using namespace xpas;


phylo_node::branch_length_type mean_branch_length(const phylo_node* root)
{
    phylo_node::branch_length_type length = 0.0;
    size_t count = 0;

    for (const auto& node : visit_subtree(root))
    {
        length += node.get_branch_length();
        count++;
    }
    return length / phylo_node::branch_length_type(count);
}

/// Inserts ghost nodes into a phylogenetic tree.
class tree_extender
{
public:
    explicit tree_extender(phylo_tree& tree)
        : _tree(tree)
          , _counter(_tree.get_node_count() + 1)
    {}
    tree_extender(const tree_extender&) = delete;
    ~tree_extender() noexcept = default;

    ghost_mapping extend()
    {
        /// add ghost nodes
        extend_subtree(_tree.get_root());

        /// reindex the tree
        _tree.index();

        /// return the mapping Extended Node ID -> Original Postorder ID
        return _mapping;
    }

private:
    void extend_subtree(phylo_node* node)
    {
        /// iterate over the copy, because we need to remove and insert elements
        /// during iteration
        auto children = node->get_children();
        for (auto child : children)
        {
            extend_subtree(child);
        }

        /// if not root
        if (node->get_parent())
        {
            auto parent = node->get_parent();

            auto old_branch_length = node->get_branch_length();

            /// The mean branch length in the subtree of the node
            auto mean_length = mean_branch_length(node);

            const auto x0_name = std::to_string(_counter++) + "_X0";
            auto x0 = new phylo_node(x0_name,  old_branch_length / 2.0, parent);
            parent->remove_child(node);
            parent->add_child(x0);

            const auto x1_name = std::to_string(_counter++) + "_X1";
            auto x1 = new phylo_node(x1_name, mean_length + old_branch_length / 2.0, x0);
            x0->add_child(x1);
            x0->add_child(node);
            node->set_branch_length(old_branch_length / 2.0);

            auto x2 = new phylo_node(std::to_string(_counter++) + "_X2",  0.01, x1);
            auto x3 = new phylo_node(std::to_string(_counter++) + "_X3",  0.01, x1);
            x1->add_child(x2);
            x1->add_child(x3);

            /// Map the new ghost node IDs to the original post-order ID.
            /// This is needed to group x0 and x1 together,
            /// and for output purposes.
            _mapping[x0_name] = node->get_postorder_id();
            _mapping[x1_name] = node->get_postorder_id();
        }
    }

    phylo_tree& _tree;
    size_t _counter;

    ghost_mapping _mapping;
};

ghost_mapping extend_tree(phylo_tree& tree)
{
    tree_extender extender(tree);
    return extender.extend();
}


std::tuple<phylo_tree, phylo_tree, ghost_mapping> xpas::preprocess_tree(const std::string& filename, bool force_root)
{
    /// load original tree
    auto tree = xpas::io::load_newick(filename);

    if (!tree.is_rooted())
    {
        if (!force_root)
        {
            throw std::runtime_error("This reference tree is not rooted."
                                     "Please provide a rooted tree or re-root it by adding --force-root. "
                                     "The trifurcation described in the newick file will be used to root the tree."
                                     "WARNING! This may impact placement accuracy.");
        }
        else
        {
            reroot_tree(tree);
        }
    }

    /// inject ghost nodes
    auto mapping = extend_tree(tree);

    auto original_tree = xpas::io::load_newick(filename);

    return std::make_tuple(std::move(original_tree), std::move(tree), mapping);
}

void xpas::reroot_tree(phylo_tree& tree)
{
    auto root = tree.get_root();

    /// If the root has 3 children
    auto children = root->get_children();
    if (children.size() > 2)
    {
        auto a = children[0];

        /// change it from (a, b, c); to ((b, c), a)added_root;
        auto new_node = new phylo_node("added_root", 0.0, nullptr);
        new_node->add_child(root);
        new_node->add_child(a);
        root->remove_child(a);

        tree.set_root(new_node);
        tree.index();
    }
}