#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main()

#include <boost/filesystem.hpp>
#include <catch2/catch.hpp>
#include <xpas/phylo_kmer_db.h>
#include <xpas/serialization.h>
#include <xpas/kmer_iterator.h>
#include <xpas/newick.h>
#include <xpas/phylo_tree.h>
#include <utils/io/fasta.h>

namespace fs = boost::filesystem;

auto create_test_map()
{
    std::unordered_map<xpas::phylo_kmer::key_type,
        std::unordered_map<xpas::phylo_kmer::branch_type, xpas::phylo_kmer::score_type>> values =
        {
            {
                0, { { 0, 0.00f } }
            },
            {
                1, { { 0, 0.10f }, { 1, 0.11f } }
            },
            {
                2, { { 0, 0.20f }, { 1, 0.21f }, { 2, 0.22f } }
            },
            {
                3, { { 1, 0.31f }, { 2, 0.32f } }
            },
            {
                4, { { 2, 0.42f } }
            }
        };
    return values;
}

template<typename MapType>
xpas::phylo_kmer_db create_db_from_map(const MapType& values, size_t kmer_size, xpas::phylo_kmer::score_type omega)
{
    xpas::phylo_kmer_db db { kmer_size, omega, xpas::seq_type::name, "" };
    for (const auto& [key, entries] : values)
    {
        for (const auto& [branch, score] : entries)
        {
            db.unsafe_insert(key, { branch, score });
        }
    }
    return db;
}

TEST_CASE("Database size", "[database]")
{
    {
        const auto values = create_test_map();
        const auto db = create_db_from_map(values, 3, 1.0);
        REQUIRE(db.size() == values.size());
    }

    {
        const xpas::phylo_kmer_db db { 3, 1.0, xpas::seq_type::name, "" };
        REQUIRE(db.size() == 0);
    }
}

TEST_CASE("K-mer size and omega", "[database]")
{
    const size_t kmer_size = 5;
    const xpas::phylo_kmer::score_type omega = 1.0;
    const xpas::phylo_kmer_db db { kmer_size, omega, xpas::seq_type::name, "" };

    REQUIRE(db.kmer_size() == kmer_size);
}

template<typename MapType>
void compare_db(const MapType& values, const xpas::phylo_kmer_db& db)
{
    for (const auto& [key, entries] : values)
    {
        auto db_entries = db.search(key);
        REQUIRE((bool)db_entries);
        for (const auto&[branch, score] : *db_entries)
        {
            REQUIRE(entries.find(branch) != entries.end());
            REQUIRE(entries.find(branch)->second == Approx(score));
        }
    }
}

TEST_CASE("Database search", "[database]")
{
    const auto filename = fs::unique_path().string();
    const auto values = create_test_map();
    const auto db = create_db_from_map(values, 3, 1.0 );
    compare_db(values, db);
}


TEST_CASE("(De-)serialization", "[database]")
{
    const auto filename = fs::unique_path().string();
    const auto values = create_test_map();
    const size_t kmer_size = 3;
    const xpas::phylo_kmer::score_type omega = 1.0;

    {
        const auto db = create_db_from_map(values, kmer_size, omega);
        xpas::save(db, filename);
    }

    {
        const auto db = xpas::load(filename);
        REQUIRE(db.size() == values.size());
        REQUIRE(db.kmer_size() == kmer_size);
        REQUIRE(db.omega() == omega);
    }
}

/// Iterate over all the combinations with repetition
/// http://shoaib-ahmed.com/2018/for-each-combination-with-repetetion-c++/
template<typename V, typename Callable>
void for_each_combination(V &v, size_t gp_sz, Callable f) {
    V gp(gp_sz);
    auto total_n = std::pow(v.size(), gp.size());
    for (auto i = 0; i < total_n; ++i) {
        auto n = i;
        for (auto j = 0ul; j < gp.size(); ++j) {
            gp[gp.size() - j - 1] = v[n % v.size()];
            n /= v.size();
        }
        f(gp);
    }
}

TEST_CASE("Encoding and decoding k-mers", "[kmers]")
{
    auto alphabet = std::vector<char>{ 'A', 'C', 'G', 'T', '-', 'N' };
    const size_t kmer_size = 3;
    size_t count = 0;
    for_each_combination(alphabet, kmer_size,
                         [&count](std::vector<char>& bases) {
                             const auto kmer = std::string{ bases.begin(), bases.end() };
                             if (const auto key = xpas::encode_kmer<xpas::no_ambiguity_policy>(kmer); key)
                             {
                                 REQUIRE(kmer == xpas::decode_kmer(*key, kmer.size()));
                                 REQUIRE(*key == count);
                                 ++count;
                             }
                         });
}

TEST_CASE("xpas::to_kmers iteration", "[kmers]")
{
    /// Simple iteration
    const auto long_read = std::string{ "--TTTAT-AAATGNNNN-CAAAN.NNTTTT---" };
    const size_t kmer_size = 4;

    size_t count = 0;
    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>(long_read, kmer_size))
    {
        REQUIRE(xpas::encode_kmer<xpas::no_ambiguity_policy>(kmer) == code);
        ++count;
    }
    REQUIRE(count == 6);
}

TEST_CASE("xpas::to_kmers empty set", "[kmers]")
{
    const size_t kmer_size = 3;

    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>("", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>("-", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>("---", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>("----", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>("NNN", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>("AAN", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>("ANA", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>("ANA", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>("NAA", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>("-AA", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>("A-A", kmer_size))
    {
        REQUIRE_FALSE(true);
    }

    for (const auto& [kmer, code] : xpas::to_kmers<xpas::no_ambiguity_policy>("AA-", kmer_size))
    {
        REQUIRE_FALSE(true);
    }
}

TEST_CASE("xpas::io::parse_newick", "[tree]")
{
    std::string newick = "((A:0.05,B:0.1):0.15,(C:0.2,D:0.25):0.3):0.35;";
    /// labels in the DFS post-order
    const std::vector<std::string> labels = { "A", "B", "", "C", "D", "", ""};
    /// branch lengths in the DFS post-order
    const std::vector<double> lengths = { 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35 };

    auto tree = xpas::io::parse_newick(newick);

    /// iterate over nodes and check if the labels and branch lengths are correct
    size_t i = 0;
    for (auto& node : tree)
    {
        REQUIRE(node.get_label() == labels[i]);
        REQUIRE(Approx(node.get_branch_length()) == lengths[i]);

        // ensure we can not create trees from non-root nodes
        if (node.get_parent())
        {
            REQUIRE_THROWS(xpas::phylo_tree{ &node });
        }

        ++i;
    }

    REQUIRE(tree.get_node_count() == 7);
}

TEST_CASE("xpas::impl::next_by_postorder", "[tree]")
{
    std::string newick = "((A:0.05,B:0.1):0.15,(C:0.2,D:0.25):0.3):0.35;";
    const auto tree = xpas::io::parse_newick(newick);

    for (const auto& node : tree)
    {
        auto next = xpas::impl::next_by_postorder(&node);

        /// next can be null only for the root node
        if (!next)
        {
            REQUIRE(node.get_parent() == nullptr);
        }
        else
        {
            REQUIRE(node.get_postorder_id() + 1 == next->get_postorder_id());
        }
    }
}

TEST_CASE("xpas::phylo_tree::get_by_postorder_id", "[tree]")
{
    std::string newick = "((A:0.05,B:0.1):0.15,(C:0.2,D:0.25):0.3):0.35;";

    const auto tree = xpas::io::parse_newick(newick);

    int iteration_count = 0;
    /// test phylo_tree::get_post_order
    for (const auto& node : tree)
    {
        // make sure this node has got the right post-order id
        REQUIRE(node.get_postorder_id() == iteration_count);

        /// find the same node in the tree by their post-order id
        const auto postorder_id = node.get_postorder_id();
        const auto found = tree.get_by_postorder_id(postorder_id);

        // make sure we found it
        REQUIRE(found);

        const auto& node_found = *found;
        // make sure it's not a null pointer
        REQUIRE(node_found);

        // compare fields
        REQUIRE(node_found->get_label() == node.get_label());
        REQUIRE(node_found->get_postorder_id() == node.get_postorder_id());
        REQUIRE(node_found->get_preorder_id() == node.get_preorder_id());
        REQUIRE(node_found->get_children() == node.get_children());
        REQUIRE(node_found->get_branch_length() == node.get_branch_length());

        ++iteration_count;
    }
}

TEST_CASE("xpas::visit_tree", "[tree]")
{
    std::string newick = "((A:0.05,B:0.1):0.15,(C:0.2,D:0.25):0.3):0.35;";
    auto tree = xpas::io::parse_newick(newick);

    std::vector<double> total_lengths = { 0.05, 0.1, 0.3, 0.2, 0.25, 0.75, 1.4 };
    size_t i = 0;

    /// can't start visiting nullptr
    REQUIRE_THROWS(xpas::visit_subtree(nullptr));

    /// Here we also test non-const iteration
    for (auto& node : tree)
    {
        /// Test if we can start visiting the subtree
        REQUIRE_NOTHROW(xpas::visit_subtree(&node));

        /// run DFS from a node, calculating the total subtree branch length
        double total_length = 0.0;
        for (const auto& subtree_node : xpas::visit_subtree(&node))
        {
            total_length += subtree_node.get_branch_length();
        }

        REQUIRE(Approx(total_length) == total_lengths[i]);

        ++i;
    }
}

TEST_CASE("xpas::io::read_fasta", "[utils]")
{
    const auto fasta = std::vector<xpas::io::fasta>{
        {"1", "AAAAA"},
        {"2", "AAAAC"},
        {"3", "AAAAG"},
        {"4", "AAAAT"},
        {"5", "AAACA"},
        {"6", "AAACC"},
        {"7", "AAACG"},
        {"8", "AAACT"},
        {"9", "AAAGA"},
        {"10", "AAAGC"},
        {"11", "AAAGG"},
        {"12", "AAAGT"},
    };

    /// write sequences to a temporary file
    const auto filename = fs::unique_path().string();
    std::ofstream out(filename);
    for (const auto& seq : fasta)
    {
        out << ">" << seq.header() << std::endl << seq.sequence() << std::endl;
    }

    /// read with default batch size
    size_t i = 0;
    for (const auto& seq : xpas::io::read_fasta(filename))
    {
        REQUIRE(seq == fasta[i]);
        ++i;
    }

    /// batch_size < fasta.size()
    /// fasta.size() % batch_size == 0
    size_t batch_size = 4;
    i = 0;
    for (const auto& seq : xpas::io::read_fasta(filename, batch_size))
    {
        REQUIRE(seq == fasta[i]);
        ++i;
    }

    /// batch_size < fasta.size()
    /// fasta.size() % batch_size != 0
    batch_size = 5;
    i = 0;
    for (const auto& seq : xpas::io::read_fasta(filename, batch_size))
    {
        REQUIRE(seq == fasta[i]);
        ++i;
    }

    /// batch_size == fasta.size()
    batch_size = fasta.size();
    i = 0;
    for (const auto& seq : xpas::io::read_fasta(filename, batch_size))
    {
        REQUIRE(seq == fasta[i]);
        ++i;
    }

    /// batch_size > fasta.size()
    batch_size = 2 * fasta.size();
    i = 0;
    for (const auto& seq : xpas::io::read_fasta(filename, batch_size))
    {
        REQUIRE(seq == fasta[i]);
        ++i;
    }
}