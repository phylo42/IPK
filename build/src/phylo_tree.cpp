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
    _init_tree();
}

phylo_tree::~phylo_tree() noexcept
{
    delete _root;
}

phylo_tree::const_iterator xpas::phylo_tree::begin() const noexcept
{
    return visit_subtree<true>(_root).begin();
}

phylo_tree::const_iterator xpas::phylo_tree::end() const noexcept
{
    return visit_subtree<true>(_root).end();
}

phylo_tree::iterator xpas::phylo_tree::begin() noexcept
{
    return visit_subtree<false>(_root).begin();
}

phylo_tree::iterator xpas::phylo_tree::end() noexcept
{
    return visit_subtree<false>(_root).end();
}

size_t phylo_tree::get_node_count() const noexcept
{
    return _node_count;
}

phylo_tree::value_pointer phylo_tree::get_root() const noexcept
{
    return _root;
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

void phylo_tree::_init_tree()
{
    if (_root->get_parent())
    {
        throw std::invalid_argument{ "Can not create a tree from non-root node: "
                                     "the parent of the root must be nullptr."};
    }

    _node_count = 0;

    auto it = visit_subtree<false>(_root).begin();
    const auto end = visit_subtree<false>(_root).end();
    for (; it != end; ++it)
    {
        const phylo_node* node{ it };
        _preorder_id_node_mapping[node->get_preorder_id()] = node;
        _postorder_id_node_mapping[node->get_postorder_id()] = node;

        ++_node_count;
    }
}

namespace xpas
{

    phylo_node::branch_length_type mean_branch_length(const phylo_node* root)
    {
        phylo_node::branch_length_type length = 0.0;
        size_t count = 0;

        for (const auto& node : visit_subtree<true>(root))
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

        void extend()
        {
            /// add ghost nodes
            extend_subtree(_tree.get_root());

            /// traverse the tree to reinitialize postorder/preorder node ids
            _tree._init_tree();
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

                auto x0 = new phylo_node(std::to_string(_counter++) + "_X0",  old_branch_length / 2.0, parent);
                parent->remove_child(node);
                parent->add_child(x0);

                auto x1 = new phylo_node(std::to_string(_counter++) + "_X1",  mean_length + old_branch_length / 2.0, x0);
                x0->add_child(x1);
                x0->add_child(node);
                node->_parent = x0;
                node->set_branch_length(old_branch_length / 2.0);

                auto x2 = new phylo_node(std::to_string(_counter++) + "_X2",  0.01, x1);
                auto x3 = new phylo_node(std::to_string(_counter++) + "_X3",  0.01, x1);
                x1->add_child(x2);
                x1->add_child(x3);
            }
        }

        phylo_tree& _tree;
        size_t _counter;

        //extended_node_mapping _extended_mapping;
    };

    void extend_tree(phylo_tree& tree)
    {
        tree_extender extender(tree);
        extender.extend();
    }
}

void save_tree(const phylo_tree& tree, const std::string& filename)
{
    std::ofstream out(filename);
    out << xpas::io::to_newick(tree);
}

phylo_tree xpas::preprocess_tree(const std::string& working_dir, const std::string& filename)
{
    /// load original tree
    auto tree = xpas::io::load_newick(filename);

    /// root if necessary
    /// ...

    /// inject ghost nodes
    extend_tree(tree);

    /// save extended tree
    fs::path extended_tree_path = fs::path(working_dir) / "extended_trees" / filename;
    save_tree(tree, extended_tree_path.string());

    return tree;
}