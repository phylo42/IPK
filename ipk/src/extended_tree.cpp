#include "extended_tree.h"
#include <i2l/newick.h>

using namespace i2l;
using namespace ipk;

phylo_node::branch_length_type total_branch_length(const phylo_node* root)
{
    if (root->is_leaf())
    {
        return 0.0;
    }
    else
    {
        phylo_node::branch_length_type length = 0.0;
        for (const auto& node : visit_subtree(root))
        {
            if (node.is_leaf())
            {
                length += node.get_branch_length();
            }
            else
            {
                length += node.get_num_leaves() * node.get_branch_length();
            }
        }
        /// exclude the branch that leads to the root, because it is not
        /// in the subtree
        length -= root->get_num_leaves() * root->get_branch_length();
        return length;
    }
}

/// Calculates the branch lengths from X0, X1 ghost nodes to their parents
std::pair<phylo_node::branch_length_type, phylo_node::branch_length_type> calc_ghost_branch_lengths(const phylo_node* node)
{
    auto old_branch_length = node->get_branch_length();

    /* The notation for different branch lengths is the following.
     * By expanding the original tree, we create new ghost nodes x0, x1:
     *
     *   parent                parent
     *      |                     |
     *      |                     |
     *      |       ====>         x0---------+
     *      |                     |          |
     *      |                     |          |
     *    node                   node        |
     *     / \                   / \         x1
     *    /   \                 /   \
     *  (subtree)             (subtree)
     *
     *
     *  This function returns the lengths of branches (x0 -> parent), (x1 -> x0)
     *
     */
    phylo_node::branch_length_type x0_branch_length = old_branch_length / 2.0;
    phylo_node::branch_length_type x1_branch_length;
    const auto residual_bl = old_branch_length - x0_branch_length;
    if (node->is_leaf())
    {
        x1_branch_length = residual_bl;
    }
    else
    {
        /// The mean branch length in the subtree of X0 can be recalculated from the
        /// mean branch length of the subtree of Node
        const auto total_length = total_branch_length(node);
        x1_branch_length = (total_length + residual_bl*node->get_num_leaves()) / node->get_num_leaves();
    }

    return {x0_branch_length, x1_branch_length };
}

/// Inserts ghost nodes into a phylogenetic tree.
class tree_extender
{
public:
    explicit tree_extender(const phylo_tree& original_tree)
        : _original_tree(original_tree)
        , _counter(original_tree.get_node_count() + 1)
    {}
    tree_extender(const tree_extender&) = delete;
    ~tree_extender() noexcept = default;

    std::pair<phylo_tree, ghost_mapping> extend()
    {
        /// copy the tree.
        auto extended_tree = _original_tree.copy();

        /// add ghost nodes
        extend_subtree(extended_tree.get_root());

        /// reindex the tree
        extended_tree.index();

        /// return the tree and
        ///        the mapping Extended Node ID -> Original Postorder ID
        return { std::move(extended_tree), std::move(_mapping) };
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

            /// Find the corresponding node in the original tree
            /// We can use the postorder IDs here because despite that the extended tree
            /// has got more nodes, it is not indexed yet, and all IDs are old
            const auto original_node = _original_tree.get_by_postorder_id(node->get_postorder_id());

            /// Use it to calculate branch length. It is easier to do in the original tree
            const auto& [x0_length, x1_length] = calc_ghost_branch_lengths(*original_node);

            const auto x0_name = std::to_string(_counter++) + "_X0";
            auto x0 = new phylo_node(x0_name, x0_length, parent);
            parent->remove_child(node);
            parent->add_child(x0);

            const auto x1_name = std::to_string(_counter++) + "_X1";

            auto x1 = new phylo_node(x1_name, x1_length, x0);
            x0->add_child(x1);
            x0->add_child(node);
            const auto old_branch_length = node->get_branch_length();
            node->set_branch_length(old_branch_length - x0_length);

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

    const phylo_tree& _original_tree;
    size_t _counter;
    ghost_mapping _mapping;
};

std::pair<phylo_tree, ghost_mapping> extend_tree(const phylo_tree& tree)
{
    tree_extender extender(tree);
    return extender.extend();
}


std::tuple<phylo_tree, phylo_tree, ghost_mapping> ipk::preprocess_tree(const std::string& filename, bool use_unrooted)
{
    /// load original tree
    auto tree = i2l::io::load_newick(filename);

    if (!tree.is_rooted())
    {
        if (!use_unrooted)
        {
            throw std::runtime_error("This reference tree is not rooted."
                                     "Please provide a rooted tree or provide --use-unrooted."
                                     "WARNING! This may impact placement accuracy.");
        }
    }

    /// inject ghost nodes
    auto [extended_tree, mapping] = extend_tree(tree);

    auto original_tree = i2l::io::load_newick(filename);
    return std::make_tuple(std::move(original_tree), std::move(extended_tree), mapping);
}

void ipk::reroot_tree(phylo_tree& tree)
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