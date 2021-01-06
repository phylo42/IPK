#include <xpas/phylo_kmer.h>
#include <phylo_tree.h>
#include <newick.h>
#include <boost/filesystem.hpp>
#include <algorithm>

namespace fs = boost::filesystem;
using namespace xpas;
using namespace xpas::impl;
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

    _preorder_id_node_mapping = std::move(other._preorder_id_node_mapping);
    _postorder_id_node_mapping = std::move(other._postorder_id_node_mapping);
}

phylo_tree::~phylo_tree() noexcept
{
    delete _root;
}

phylo_tree::const_iterator xpas::phylo_tree::begin() const noexcept
{
    return visit_subtree(_root).begin();
}

phylo_tree::const_iterator xpas::phylo_tree::end() const noexcept
{
    return visit_subtree(_root).end();
}

phylo_tree::iterator xpas::phylo_tree::begin() noexcept
{
    return visit_subtree<postorder_tree_iterator<false>>(_root).begin();
}

phylo_tree::iterator xpas::phylo_tree::end() noexcept
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

void phylo_tree::index()
{
    if (_root->get_parent())
    {
        throw std::invalid_argument{ "Can not create a tree from non-root node: "
                                     "the parent of the root must be nullptr."};
    }


    /// Set postorder id for all nodes
    phylo_node::id_type postorder_id = 0;
    for (auto& node : xpas::visit_subtree<postorder_tree_iterator<false>>(_root))
    {
        node._postorder_id = postorder_id;
        ++postorder_id;
    }

    /// Set preorder id for all nodes
    phylo_node::id_type preorder_id = 0;

    //auto it = preorder_tree_iterator<false>(_root);
    //const auto end = preorder_tree_iterator<false>(nullptr);
    //while (it != end)
    for (auto& node : xpas::visit_subtree<preorder_tree_iterator<false>>(_root))
    {
        node._preorder_id = preorder_id;
        ++preorder_id;
    }

    /// Make maps for fast search by postorder/preorder id
    _preorder_id_node_mapping.clear();
    _postorder_id_node_mapping.clear();
    _node_count = 0;

    for (const auto& node : xpas::visit_subtree(_root))
    {
        _preorder_id_node_mapping[node.get_preorder_id()] = &node;
        _postorder_id_node_mapping[node.get_postorder_id()] = &node;

        ++_node_count;
    }
}

optional<const phylo_node*> phylo_tree::get_by_preorder_id(phylo_node::id_type preorder_id) const noexcept
{
    if (const auto it = _preorder_id_node_mapping.find(preorder_id); it != _postorder_id_node_mapping.end())
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

namespace xpas
{

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
                node->_parent = x0;
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
}

void xpas::save_tree(const phylo_tree& tree, const std::string& filename)
{
    std::ofstream out(filename);
    out << xpas::io::to_newick(tree);
}


std::tuple<phylo_tree, phylo_tree, ghost_mapping> xpas::preprocess_tree(const std::string& filename)
{
    /// load original tree
    auto tree = xpas::io::load_newick(filename);

    /// root if necessary
    /// ...

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