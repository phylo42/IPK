#include <core/phylo_tree.h>
#include <core/phylo_kmer.h>

using namespace core;
using namespace core::impl;
using std::vector;
using std::string;
using std::move;
using std::begin, std::end;

phylo_node::phylo_node()
{
    _clean();
}

phylo_node::phylo_node(id_type postorder_id, const string& label, float branch_length,
    const vector<phylo_node*>& children, phylo_node* parent)
    : _postorder_id{ postorder_id }
    , _label{ label }
    , _branch_length{ branch_length }
    , _children{ children }
    , _parent{ parent }
{
}

phylo_node::~phylo_node() noexcept
{
    for (auto child : _children)
    {
        delete child;
    }
}

bool phylo_node::operator==(const phylo_node& rhs) const noexcept
{
    return (_postorder_id == rhs._postorder_id) && (_label == rhs._label);
}

bool phylo_node::operator!=(const phylo_node& rhs) const noexcept
{
    return !operator==(rhs);
}

std::string phylo_node::get_label() const
{
    return _label;
}

phylo_node* phylo_node::get_parent() const
{
    return _parent;
}

float phylo_node::get_branch_length() const
{
    return _branch_length;
}

std::vector<phylo_node*> phylo_node::get_children() const
{
    return _children;
}

void phylo_node::_clean()
{
    _postorder_id = -1;
    _label = "";
    _branch_length = 0;
    _children.clear();
    _parent = nullptr;
}

void phylo_node::_add_children(phylo_node* node)
{
    _children.push_back(node);
}

core::phylo_tree::phylo_tree(phylo_node* root, size_t node_count) noexcept
    : _root{ root }, _node_count{ node_count }
{}

core::phylo_tree::~phylo_tree() noexcept
{
    delete _root;
}

phylo_node* core::impl::get_leftmost_leaf(phylo_node* root)
{
    while (!root->get_children().empty())
    {
        root = root->get_children()[0];
    }
    return root;
}

core::phylo_tree::const_iterator core::phylo_tree::begin() const
{
    return postorder_tree_iterator<true>{ core::impl::get_leftmost_leaf(_root) };
}

core::phylo_tree::const_iterator core::phylo_tree::end() const
{
    return postorder_tree_iterator<true>(nullptr);
}

size_t core::phylo_tree::get_node_count() const
{
    return _node_count;
}
