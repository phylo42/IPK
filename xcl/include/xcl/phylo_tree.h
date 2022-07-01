#ifndef XPAS_PHYLO_TREE_H
#define XPAS_PHYLO_TREE_H

#include <xcl/optional.h>
#include <xcl/phylo_kmer.h>
#include <xcl/phylo_node.h>

namespace xpas {

    namespace impl
    {
        /// \brief Finds the leftmost leaf of a subtree.
        /// \details Used to start a depth-first search.
        /// The Const template parameter determines whether the visited nodes are returned via const
        /// or non-const references and pointers.
        template<bool Const = true,
            typename pointer = std::conditional_t<Const, const phylo_node*, phylo_node*>>
        pointer get_leftmost_leaf(pointer root) noexcept
        {
            while (root && !root->get_children().empty())
            {
                root = root->get_children()[0];
            }
            return root;
        }


        /// \brief Find the index of this node in the parent's array of children
        /// \details The Const template parameter determines whether the visited nodes are returned via const
        /// or non-const references and pointers.
        template<bool Const = true,
            typename pointer = std::conditional_t<Const, const phylo_node*, phylo_node*>>
        optional<phylo_node::id_type> _id_in_parent(pointer node)
        {
            if (node && node->get_parent() != nullptr)
            {
                /// WARNING:
                /// Here we perform a linear search to look for an index
                /// in the parent's children list.
                const auto& children = node->get_parent()->get_children();
                const auto it = std::find(begin(children), end(children), node);
                if (it != end(children))
                {
                    return {distance(begin(children), it)};
                }
            }
            return {nullopt};
        }

        /// \brief Return the node that would have the next post-order id in the tree,
        /// if such node exists. Return null pointer otherwise.
        template<bool Const = true,
            typename pointer = std::conditional_t<Const, const phylo_node*, phylo_node*>>
        pointer next_by_postorder(pointer node)
        {
            /// Go upside down if necessary. We need to know the index of current
            /// node in the parent->children
            pointer temp = node->get_parent();
            auto idx = _id_in_parent(node);
            while (!idx && temp)
            {
                temp = node->get_parent();
                idx = _id_in_parent(node);
            }

            /// check if we have reached the end of the tree
            if (temp == nullptr)
            {
                node = nullptr;
            }
            /// if not, and there are siblings, visit the next sibling
            else if (auto sibling_id = static_cast<size_t>(*idx + 1); sibling_id < temp->get_children().size())
            {
                node = temp->get_children()[sibling_id];
                node = get_leftmost_leaf(node);
            }
            /// otherwise visit the parent
            else
            {
                node = temp;
            }
            return node;
        }

        /// \brief Return the node that would have the next pre-order id in the tree,
        /// if such node exists. Return null pointer otherwise.
        template<bool Const = true,
            typename pointer = std::conditional_t<Const, const phylo_node*, phylo_node*>>
        pointer next_by_preorder(pointer node)
        {
            /// if there are children, go to the left subtree
            if (node->get_children().size())
            {
                return node->get_children()[0];
            }
            else
            {
                pointer parent = node->get_parent();
                auto idx = _id_in_parent(node);

                /// go up until we find an unvisited sibling
                auto sibling_id = static_cast<size_t>(*idx + 1);
                while (parent && sibling_id >= parent->get_children().size())
                {
                    node = parent;
                    parent = node->get_parent();
                    idx = _id_in_parent(node);
                    sibling_id = static_cast<size_t>(*idx + 1);
                }

                /// if found a sibling, visit it
                if (parent)
                {
                    return parent->get_children()[sibling_id];
                }
                /// otherwise it's the end of the tree
                else
                {
                    return nullptr;
                }
            }
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief A forward access non-const iterator for phylo_node objects.
    /// \details Performs a pre-order depth-first search among a subtree of an input phylo_node.
    /// The Const template parameter determines whether the visited nodes are returned via const
    /// or non-const references and pointers.
    template<bool Const=true>
    class preorder_tree_iterator
    {
    public:
        static constexpr bool value = Const;

        using iterator_category = std::forward_iterator_tag;

        /// deduce const qualifier from bool Const parameter
        using reference = typename std::conditional_t<Const, const phylo_node&, phylo_node&>;
        using pointer = typename std::conditional_t<Const, const phylo_node*, phylo_node*>;

        /// How to start and stop the iterator
        static constexpr auto begin = [](pointer root) { return root; };
        static constexpr auto end = [](pointer) { return nullptr; };
    public:

        preorder_tree_iterator() noexcept
            : preorder_tree_iterator{ nullptr }
        {}

        explicit preorder_tree_iterator(pointer node) noexcept
            : _current{ node }
        {}

        preorder_tree_iterator& operator=(const preorder_tree_iterator& rhs) noexcept
        {
            if (*this != rhs)
            {
                _current = rhs._current;
            }
            return *this;
        }

        explicit operator pointer() const noexcept
        {
            return _current;
        }

        bool operator==(const preorder_tree_iterator& rhs) const noexcept
        {
            return _current == rhs._current;
        }

        bool operator!=(const preorder_tree_iterator& rhs) const noexcept
        {
            return !(*this == rhs);
        }

        preorder_tree_iterator& operator++()
        {
            _current = impl::next_by_preorder(_current);
            return *this;
        }

        reference operator*() const
        {
            return *_current;
        }

        pointer operator->() const noexcept
        {
            return _current;
        }

    private:
        pointer _current;
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief A forward access non-const iterator for phylo_node objects.
    /// \details Performs a postorder depth-first search among a subtree of an input phylo_node.
    /// The Const template parameter determines whether the visited nodes are returned via const
    /// or non-const references and pointers.
    template<bool Const=true>
    class postorder_tree_iterator
    {
    public:
        static constexpr bool value = Const;

        using iterator_category = std::forward_iterator_tag;

        /// deduce const qualifier from bool Const parameter
        using reference = typename std::conditional_t<Const, const phylo_node&, phylo_node&>;
        using pointer = typename std::conditional_t<Const, const phylo_node*, phylo_node*>;

        /// How to start and stop the iterator
        static constexpr auto begin = impl::get_leftmost_leaf<Const>;
        static constexpr auto end = impl::next_by_postorder<Const>;
    public:

        postorder_tree_iterator() noexcept
            : postorder_tree_iterator{ nullptr }
        {}

        explicit postorder_tree_iterator(pointer node) noexcept
                : _current{ node }
        {}

        postorder_tree_iterator& operator=(const postorder_tree_iterator& rhs) noexcept
        {
            if (*this != rhs)
            {
                _current = rhs._current;
            }
            return *this;
        }

        explicit operator pointer() const noexcept
        {
            return _current;
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
            _current = impl::next_by_postorder(_current);
            return *this;
        }

        reference operator*() const
        {
            return *_current;
        }

        pointer operator->() const noexcept
        {
            return _current;
        }

    private:
        pointer _current;
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief A wrapper for postorder_tree_iterator. Visits subtree nodes in DFS post-order.
    /// \details The Const template paramter has the same meaning as in impl::postorder_tree_iterator.
    template<class Iterator=postorder_tree_iterator<true>,
             bool Const=Iterator::value>
    class visit_subtree
    {
    public:
        using iterator = Iterator;
        using value_type = phylo_node;

        /// deduce const qualifier from bool Const parameter
        using reference = typename std::conditional_t<Const, const phylo_node&, phylo_node&>;
        using pointer = typename std::conditional_t<Const, const phylo_node*, phylo_node*>;

        explicit visit_subtree(pointer root)
            : _root { root }
        {
            /// check if it is a null pointer, or not a
            if (!root)
            {
                throw std::runtime_error("xpas::visit_subtree: root can not be null");
            }
        }

        [[nodiscard]]
        iterator begin() const noexcept
        {
            return iterator{ Iterator::begin(_root) };
        }

        [[nodiscard]]
        iterator end() const noexcept
        {
            return iterator{ Iterator::end(_root) };
        }

    private:
        pointer _root;
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief A phylogenetic tree class
    /// \defails phylo_tree is only constructable by the rappas::io::load_newick function.
    /// Non-copyable. Phylo-nodes are not modifiable.
    /// \sa xpas::phylo_node, xpas::io::load_newick
    class phylo_tree
    {
        friend class tree_extender;
    public:
        /// Member types
        /// const iterator type. Performs a post-order depth-first search
        using const_iterator = postorder_tree_iterator<true>;
        /// non-const iterator type.  Performs a post-order depth-first search
        using iterator = postorder_tree_iterator<false>;
        using value_pointer = phylo_node*;

        /// Ctors, dtor and operator=
        phylo_tree(value_pointer root);
        phylo_tree(phylo_tree&&) noexcept;
        phylo_tree(const phylo_tree&) = delete;
        phylo_tree& operator=(const phylo_tree&) = delete;
        phylo_tree& operator=(phylo_tree&&) = delete;
        ~phylo_tree() noexcept;


        /// Iterators
        /// Returns a const iterator to the beginning
        const_iterator begin() const noexcept;
        /// Returns a const iterator to the end
        const_iterator end() const noexcept;
        /// Returns a non-const iterator to the beginning
        iterator begin() noexcept;
        /// Returns a non-const iterator to the end
        iterator end() noexcept;


        /// Access
        /// \brief Returns the number of nodes in a tree.
        size_t get_node_count() const noexcept;
        /// \brief Returns a pointer to the root node.
        value_pointer get_root() const noexcept;
        /// \brief Sets the new root
        void set_root(value_pointer root);
        /// \brief Returns true if the root has no more than two children
        bool is_rooted() const noexcept;

        /// \brief Runs DFS from the root to make pre- and post-order ID mappings
        /// for all nodes
        void index();

        /// \brief Returns an optional fora pointer to a phylo_node with a given preorder_id, if exists in the tree.
        /// \details This operation does not require any traversal and implemented in O(1).
        /// \sa get_by_preorder_id
        optional<const xpas::phylo_node*> get_by_preorder_id(phylo_node::id_type preorder_id) const noexcept;

        /// \brief Returns an optional fora pointer to a phylo_node with a given postorder_id, if exists in the tree.
        /// \details This operation does not require any traversal and implemented in O(1).
        /// \sa get_by_postorder_id
        optional<const xpas::phylo_node*> get_by_postorder_id(phylo_node::id_type postorder_id) const noexcept;

        /// \brief Returns an optional for a pointer to a phylo_node with a given label, if exists in the tree.
        /// \details This operation does not require any traversal and implemented in O(1).
        optional<const xpas::phylo_node*> get_by_label(const std::string& label) const noexcept;

        /// Creates a copy of the tree. We prefer to have this method and the copy constructor deleted
        /// to make sure we never copy by mistake, and always move trees.
        phylo_tree copy() const;
    private:
        void _index_preorder_id();
        void _index_postorder_id();
        void _index_labels();
        void _index_nodes();

        /// \brief A root node.
        value_pointer _root;

        /// \brief A total number of nodes in a tree.
        /// \details WARNING: This class does not double check if the tree actually have this
        /// number of nodes, it is just passed as an argument in the constructor.
        size_t _node_count;

        /// Map "phylo_node_preorder_id -> phylo_node"
        std::unordered_map<phylo_node::id_type, const xpas::phylo_node*> _preorder_id_to_node;
        /// Map "phylo_node_postorder_id -> phylo_node"
        std::unordered_map<phylo_node::id_type, const xpas::phylo_node*> _postorder_id_node_mapping;
        /// Map "node label -> phylo_node"
        std::unordered_map<std::string, const xpas::phylo_node*> _label_to_node;
    };

    void save_tree(const phylo_tree& tree, const std::string& filename);
}

#endif