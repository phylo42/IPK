#include <xcl/phylo_kmer.h>
#include <xcl/phylo_tree.h>
#include <xcl/newick.h>
#include <boost/filesystem.hpp>
#include <algorithm>

namespace fs = boost::filesystem;
using namespace xcl;
using namespace xcl::impl;
using std::vector;
using std::string;
using std::move;
using std::begin, std::end;


phylo_tree::phylo_tree(phylo_node* root)
    : _root{ root }, _node_count{ 0 }
{
    index();
}

phylo_tree::phylo_tree(phylo_tree&& other) noexcept
{
    _root = other._root;
    other._root = nullptr;

    _node_count = other._node_count;
    other._node_count = 0;

    _preorder_id_to_node = std::move(other._preorder_id_to_node);
    _postorder_id_node_mapping = std::move(other._postorder_id_node_mapping);
    _label_to_node = std::move(other._label_to_node);
}

phylo_tree::~phylo_tree() noexcept
{
    delete _root;
}

phylo_tree::const_iterator phylo_tree::begin() const noexcept
{
    return visit_subtree(_root).begin();
}

phylo_tree::const_iterator phylo_tree::end() const noexcept
{
    return visit_subtree(_root).end();
}

phylo_tree::iterator phylo_tree::begin() noexcept
{
    return visit_subtree<postorder_tree_iterator<false>>(_root).begin();
}

phylo_tree::iterator phylo_tree::end() noexcept
{
    return visit_subtree<postorder_tree_iterator<false>>(_root).end();
}

size_t phylo_tree::get_node_count() const noexcept
{
    return _node_count;
}

phylo_tree::value_pointer phylo_tree::get_root() const noexcept
{
    return _root;
}

void phylo_tree::set_root(value_pointer root)
{
    _root = root;
}

bool phylo_tree::is_rooted() const noexcept
{
    return _root && _root->get_children().size() < 3;
}

void phylo_tree::index()
{
    if (_root->get_parent())
    {
        throw std::invalid_argument{ "Can not create a tree from non-root node: "
                                     "the parent of the root must be nullptr."};
    }

    /// Recreate all search maps
    _index_preorder_id();
    _index_postorder_id();
    _index_labels();

    _index_nodes();
}

optional<const phylo_node*> phylo_tree::get_by_preorder_id(phylo_node::id_type preorder_id) const noexcept
{
    if (const auto it = _preorder_id_to_node.find(preorder_id); it != _postorder_id_node_mapping.end())
    {
        return { it->second };
    }
    else
    {
        return { nullopt };
    }
}

optional<const phylo_node*> phylo_tree::get_by_postorder_id(phylo_node::id_type postorder_id) const noexcept
{
    if (const auto it = _postorder_id_node_mapping.find(postorder_id); it != _postorder_id_node_mapping.end())
    {
        return { it->second };
    }
    else
    {
        return { nullopt };
    }
}

optional<const phylo_node*> phylo_tree::get_by_label(const std::string& label) const noexcept
{
    if (const auto it = _label_to_node.find(label); it != _label_to_node.end())
    {
        return { it->second };
    }
    else
    {
        return { nullopt };
    }
}

phylo_tree phylo_tree::copy() const
{
    auto new_root = _root->copy();
    return phylo_tree(new_root);
}

void phylo_tree::_index_preorder_id()
{
    _preorder_id_to_node.clear();

    phylo_node::id_type preorder_id = 0;
    for (auto& node : visit_subtree<preorder_tree_iterator<false>>(_root))
    {
        node._preorder_id = preorder_id;
        _preorder_id_to_node[preorder_id] = &node;
        ++preorder_id;
    }
}

void phylo_tree::_index_postorder_id()
{
    _postorder_id_node_mapping.clear();

    phylo_node::id_type postorder_id = 0;
    for (auto& node : visit_subtree<postorder_tree_iterator<false>>(_root))
    {
        node._postorder_id = postorder_id;
        _postorder_id_node_mapping[postorder_id] = &node;
        ++postorder_id;
    }
}

void phylo_tree::_index_labels()
{
    _label_to_node.clear();

    for (auto& node : visit_subtree(_root))
    {
        _label_to_node[node.get_label()] = &node;
    }
}

void phylo_tree::_index_nodes()
{
    _node_count = 0;

    for (auto& node : visit_subtree<iterator>(_root))
    {
        (void)node;
        ++_node_count;

        if (node.is_leaf())
        {
            /// By convention the leaves have no nodes in their subtrees
            node.set_num_nodes(0);

            /// By convention the number of leaves in the subtree of a leaf is one
            node.set_num_leaves(1);

            /// There is no subtree
            node.set_subtree_branch_length(0.0);
        }
        else
        {
            /// The number of the nodes in the subtree of this node
            ///   =  the number of children
            size_t total_num_nodes = node.get_children().size();
            size_t total_num_leaves = 0;

            /// The total branch length in the subtree
            phylo_node::branch_length_type subtree_branch_length = 0.0;

            for (const auto& child : node.get_children())
            {
                ///  + the total number of nodes in subtrees of children
                total_num_nodes += child->get_num_nodes();

                /// The number of leaves is the sum of leaves of all children
                total_num_leaves += child->get_num_leaves();

                subtree_branch_length += child->get_subtree_branch_length() + child->get_branch_length();
            }

            node.set_num_nodes(total_num_nodes);
            node.set_num_leaves(total_num_leaves);
            node.set_subtree_branch_length(subtree_branch_length);
        }
    }
}

namespace xcl
{
    void save_tree(const phylo_tree& tree, const std::string& filename)
    {
        std::ofstream out(filename);
        out << xcl::io::to_newick(tree);
    }
}