#include <xpas/phylo_tree.h>
#include <xpas/newick.h>
#include <iostream>
#include <boost/filesystem.hpp>

void write_test_tree(const std::string& filename)
{
    std::ofstream out(filename);
    out << "(A:0.1,B:0.2,((C:0.1,D:0.2)Y:0.1,(E:0.1,F:0.2)Z:0.2)X:0.2)W:0.0;";
}

int main()
{
    const auto filename = boost::filesystem::unique_path().string();
    write_test_tree(filename);

    /// const node iteration
    {
        const auto tree = xpas::io::load_newick(filename);
        std::cout << "Nodes: " << tree.get_node_count() << '\n';
        std::cout << xpas::io::to_newick(tree) <<  "\n\n";
    }

    /// non-const node iteration
    {
        auto tree = xpas::io::load_newick(filename);
        for (auto& node : tree)
        {
            node.set_label(node.get_label() + "+");
            node.set_branch_length(1.0 + node.get_branch_length());
        }

        std::cout << tree << std::endl;
    }

    /// visit const subtree
    {
        const auto tree = xpas::io::load_newick(filename);

        /// let's find a node by its post-order id
        const size_t postorder_id = 4;

        /// if found, visit the subtree
        if (const auto& node = tree.get_by_postorder_id(postorder_id); node)
        {
            size_t visited = 0;
            for (const auto& subtree_node : xpas::visit_subtree<true>(*node))
            {
                std::cout << subtree_node.get_label() << " : " << subtree_node.get_branch_length() << std::endl;
                ++visited;
            }

            std::cout << "Visited " << visited << " nodes." << std::endl;
        }
        else
        {
            std::cerr << "Could not find node by post-order id: " << postorder_id << std::endl;
        }
    }
}