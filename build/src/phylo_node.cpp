#include <xpas/phylo_kmer.h>
#include "phylo_node.h"

using namespace xpas;
using std::vector;
using std::string;
using std::move;
using std::begin, std::end;


phylo_node::phylo_node()
{
    clean();
}

phylo_node::phylo_node(std::string label, branch_length_type branch_length, phylo_node* parent)
    : _preorder_id(-1)
    , _postorder_id(-1)
    , _label(std::move(label))
    , _branch_length(branch_length)
    , _parent(parent)
{}

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

void phylo_node::set_parent(phylo_node* parent)
{
    _parent = parent;
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

void phylo_node::clean()
{
    _preorder_id = -1;
    _postorder_id = -1;
    _label = "";
    _branch_length = 0.0;
    _children.clear();
    _parent = nullptr;
}

void phylo_node::add_child(phylo_node* node)
{
    _children.push_back(node);
    node->set_parent(this);
}

void phylo_node::remove_child(phylo_node* node)
{
    _children.erase(std::remove(_children.begin(), _children.end(), node));
}

bool phylo_node::is_leaf() const noexcept
{
    return _children.empty();
}

bool phylo_node::is_root() const noexcept
{
    return _parent == nullptr;
}