#ifndef RAPPAS_CPP_PHYLO_TREE_H
#define RAPPAS_CPP_PHYLO_TREE_H

#include <string>
#include <vector>
#include <algorithm>
#include <meta.h>

namespace core
{
    class phylo_tree;
}

namespace rappas
{
    namespace io
    {
        class newick_parser;
        core::phylo_tree load_newick(const std::string& file_name);
    }
}

namespace core
{
    namespace impl
    {
        template<bool IsConst>
        class postorder_tree_iterator;

        /// \brief A node of a phylogenetic tree.
        class phylo_node
        {
            friend rappas::io::newick_parser;

            template<bool IsConst> friend
            class phylo_tree_iterator;

        public:
            phylo_node();

            explicit phylo_node(int id, const std::string& label, float branch_length,
                                const std::vector<phylo_node*>& children, phylo_node* parent);
            phylo_node(const phylo_node& other) = delete;
            phylo_node& operator=(const phylo_node&) = delete;
            ~phylo_node() noexcept;

            /// WARNING: this operator only checks for the id and label fields
            bool operator==(const phylo_node& rhs) const noexcept;
            bool operator!=(const phylo_node& rhs) const noexcept;

            std::string get_label() const;
            phylo_node* get_parent() const;
            float get_branch_length() const;

            std::vector<phylo_node*> get_children() const;

        private:
            /// Clean node and fill with the default values. Used in the default constructor
            void _clean();

            void _add_children(phylo_node* node);

        private:
            // TODO: test if _id convention is the same as RAPPAS' exented_tree_id one
            int _id;

            std::string _label;
            float _branch_length;
            std::vector<phylo_node*> _children;
            phylo_node* _parent;
        };

        phylo_node* get_leftmost_child(phylo_node* root);

        /// \brief A forward access (non-)const iterator for phylo_node objects. Performs a depth-first
        /// search among a subtree of input phylo_node.
        template<bool IsConst>
        class postorder_tree_iterator
        {
        public:
            using iterator_category = std::forward_iterator_tag;
            using reference = typename rappas::choose<IsConst, const phylo_node&, phylo_node&>::type;
            using pointer = typename rappas::choose<IsConst, const phylo_node*, phylo_node*>::type;

        public:
            postorder_tree_iterator()
                : postorder_tree_iterator(nullptr)
            {
            }

            explicit postorder_tree_iterator(phylo_node* node)
                : _current(node)
            {
            }

            postorder_tree_iterator(const postorder_tree_iterator& other) = default;
            ~postorder_tree_iterator() noexcept = default;

            postorder_tree_iterator& operator=(const postorder_tree_iterator& rhs)
            {
                if (*this != rhs)
                {
                    _current = rhs._current;
                }
                return *this;
            }

            bool operator==(const postorder_tree_iterator& rhs) const noexcept
            {
                return _current == rhs._current;
            }

            bool operator!=(const postorder_tree_iterator& rhs) const noexcept
            {
                return !(*this == rhs);
            }

            postorder_tree_iterator& operator++()
            {
                /// Go upside down if necessary. We need to know the index of current node in the parent->children
                phylo_node* temp = _current->_parent;
                int idx = _id_in_parent(_current);
                while (idx == -1 && temp)
                {
                    temp = _current->_parent;
                    idx = _id_in_parent(_current);
                }

                /// the end of the tree
                if (temp == nullptr)
                {
                    _current = nullptr;
                }
                /// visit the next sibling
                else if ((size_t) idx + 1 < temp->_children.size())
                {
                    _current = temp->_children[idx + 1];
                    _current = get_leftmost_child(_current);
                }
                /// visit the parent
                else
                {
                    _current = temp;
                }
                return *this;
            }

            reference operator*() const
            {
                return *_current;
            }

            pointer operator->() const
            {
                return _current;
            }

        private:
            int _id_in_parent(phylo_node* node) const
            {
                if (node->_parent != nullptr)
                {
                    const auto& children = node->_parent->_children;
                    const auto it = std::find(begin(children), end(children), node);
                    if (it != end(children))
                    {
                        return distance(begin(children), it);
                    }
                }
                return -1;
            }

        private:
            phylo_node* _current;
        };
    }

    /// \brief A phylogenetic tree class
    /// \defails phylo_tree is only constructable with a factory function. Also non-copyable and non-movable.
    /// \sa load_newick
    class phylo_tree
    {
        friend phylo_tree rappas::io::load_newick(const std::string& file_name);

    public:
        using const_iterator = impl::postorder_tree_iterator<true>;

        phylo_tree(phylo_tree&&) = delete;
        phylo_tree(const phylo_tree&) = delete;
        phylo_tree& operator=(const phylo_tree&) = delete;
        phylo_tree& operator=(phylo_tree&&) = delete;
        ~phylo_tree() noexcept;

        const_iterator begin() const;
        const_iterator end() const;

        size_t get_node_count() const;

    private:
        phylo_tree(impl::phylo_node* root, size_t node_count) noexcept;

        impl::phylo_node* _root;
        size_t _node_count;
    };
}

namespace rappas
{
    namespace io
    {
        /// \brief Reads a .newick formatted file
        core::phylo_tree load_newick(const std::string& file_name);
    }
}
#endif