#include <xpas/phylo_tree.h>
#include <xpas/phylo_kmer.h>
#include <algorithm>

using namespace xpas;
using namespace xpas::impl;
using std::vector;
using std::string;
using std::move;
using std::begin, std::end;

phylo_node::phylo_node()
{
    _clean();
}

phylo_node::~phylo_node() noexcept
{
    for (phylo_node* child : _children)
    {
        delete child;
    }
}

bool phylo_node::operator==(const phylo_node& rhs) const noexcept
{
    return (_preorder_id == rhs._preorder_id) &&
           (_postorder_id == rhs._postorder_id) &&
           (_label == rhs._label);
}

bool phylo_node::operator!=(const phylo_node& rhs) const noexcept
{
    return !operator==(rhs);
}

std::string phylo_node::get_label() const noexcept
{
    return _label;
}

void phylo_node::set_label(const std::string& label)
{
    _label = label;
}

phylo_node* phylo_node::get_parent() const noexcept
{
    return _parent;
}

phylo_node::id_type phylo_node::get_preorder_id() const noexcept
{
    return _preorder_id;
}

phylo_node::id_type phylo_node::get_postorder_id() const noexcept
{
    return _postorder_id;
}

phylo_node::branch_length_type phylo_node::get_branch_length() const noexcept
{
    return _branch_length;
}

void phylo_node::set_branch_length(branch_length_type length)
{
    _branch_length = length;
}

const std::vector<phylo_node*>& phylo_node::get_children() const
{
    return _children;
}

void phylo_node::_clean()
{
    _preorder_id = -1;
    _postorder_id = -1;
    _label = "";
    _branch_length = 0.0;
    _children.clear();
    _parent = nullptr;
}

void phylo_node::_add_children(phylo_node* node)
{
    _children.push_back(node);
}

phylo_tree::phylo_tree(phylo_node* root)
    : _root{ root }, _node_count{ 0 }
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