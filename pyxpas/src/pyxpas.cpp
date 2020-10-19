#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/numpy/dtype.hpp>
#include <xpas/phylo_kmer_db.h>
#include <xpas/serialization.h>

using namespace xpas;
namespace bp = boost::python;
namespace bn = boost::python::numpy;

/// Custom exceptions
struct AttributeError: std::exception
{
    [[nodiscard]]
    const char* what() const noexcept override { return "AttributeError exception"; }
};

struct TypeError: std::exception
{
    [[nodiscard]]
    const char* what() const noexcept override { return "TypeError exception"; }
};

/// exception translator
void translate(const std::exception& e)
{
    if (dynamic_cast<const AttributeError*>(&e))
    {
        PyErr_SetString(PyExc_AttributeError, e.what());
    }

    if (dynamic_cast<const TypeError*>(&e))
    {
        PyErr_SetString(PyExc_TypeError, e.what());
    }
}

/// Convert std::optional<search_result> to a python object
struct to_python_search_result
{
    static PyObject* convert(const std::optional<xpas::impl::search_result>& obj)
    {
        if (obj)
        {
            return bp::incref(bp::object(*obj).ptr());
        }
        else
        {
            /// We can throw an exception if nothing found, but it is better to
            /// create an empty vector of entries to pass to python side as a result.
            phylo_kmer_db::value_type empty_collection;
            xpas::impl::search_result empty_result(empty_collection.begin(), empty_collection.end());

            return bp::incref(bp::object(empty_result).ptr());
        }
    }
};

struct db_entry
{
    db_entry(xpas::phylo_kmer_db::key_type key, std::vector<xpas::pkdb_value> values)
        : _key {key}, _values{ std::move(values) } {}

    xpas::phylo_kmer_db::key_type _key;
    std::vector<xpas::pkdb_value> _values;
};

bp::tuple db_entry_iter(db_entry& entry)
{
    Py_intptr_t shape[1] = { static_cast<Py_intptr_t>(entry._values.size()) };

    const auto npuint = bn::dtype::get_builtin<xpas::phylo_kmer::branch_type>();
    auto branches = bn::zeros(1, shape, npuint);

    const auto npfloat = bn::dtype::get_builtin<float>();
    auto scores = bn::zeros(1, shape, npfloat);

    size_t i = 0;
    for (const auto& value : entry._values)
    {
        branches[i] = value.branch;
        scores[i] = value.score;
        ++i;
    }
    return bp::make_tuple(entry._key, branches, scores);
}

bn::ndarray db_as_matrix(const phylo_kmer_db& db)
{
    /// find the number of branches
    size_t num_branches = 0;
    for (const auto& [key, entries] : db)
    {
        for (const auto& [branch, score] : entries)
        {
            num_branches = std::max(num_branches, static_cast<size_t>(branch) + 1);
        }
    }

    /// find the number of k-mers
    const size_t max_key = std::pow(xpas::seq_traits::alphabet_size, db.kmer_size());

    const auto shape = bp::make_tuple(max_key, num_branches);
    const auto np_score_type = bn::dtype::get_builtin<xpas::phylo_kmer::score_type>();
    auto matrix = bn::zeros(shape, np_score_type);

    for (const auto& [key, entries] : db)
    {
        for (const auto& [branch, score] : entries)
        {
            matrix[key][branch] = score;
        }
    }

    return matrix;
}


struct to_python_db_entry
{
    static PyObject* convert(const phylo_kmer_db::storage::const_iterator::value_type& entry)
    {
        db_entry copy{ entry.first, entry.second };
        /*std::cout << entry.first;
        for (const auto& [x, y] : entry.second)
        {
            std::cout << " (" << x << " " << y << " )";
        };
        std::cout << std::endl;*/

        return bp::incref(bp::object(copy).ptr());
        //return bp::incref(bp::object(db_entry).ptr());
    }
};


BOOST_PYTHON_MODULE(pyxpas)
{
    Py_Initialize();
    bn::initialize();

    /// Exception translator
    bp::register_exception_translator<AttributeError>(&translate);
    bp::register_exception_translator<TypeError>(&translate);

    /// Python wrapping of std::optional.
    /// This is an adaption of the solution from stackoverlow.
    /// See https://stackoverflow.com/questions/36485840/wrap-boostoptional-using-boostpython
    bp::to_python_converter<std::optional<xpas::impl::search_result>, to_python_search_result>();

    /// a pair { branch -> score }
    bp::class_<pkdb_value>("PkDbValue", bp::init<phylo_kmer::branch_type, phylo_kmer::score_type>())
        .def_readonly("branch", &pkdb_value::branch)
        .def_readonly("score", &pkdb_value::score)
    ;

    /// a search result: collection of pkdb_value
    bp::class_<xpas::impl::search_result>("PkDbSearchResult")
        .def("__iter__", bp::range(&xpas::impl::search_result::begin, &xpas::impl::search_result::end))
    ;

    /*bp::class_<phylo_kmer_db::storage::const_iterator::value_type>( "PkDbEntry", )
        //.def("key", &phylo_kmer_db::storage::const_iterator::value_type::first)
        //.def("value", &phylo_kmer_db::storage::const_iterator::value_type::second)
    ;*/
    bp::class_<db_entry>("PkDbEntry", bp::init<xpas::phylo_kmer_db::key_type, std::vector<xpas::pkdb_value>>())
        .def("as_row", &db_entry_iter)
    ;
    bp::to_python_converter<phylo_kmer_db::storage::const_iterator::value_type, to_python_db_entry>();

    /// phylo k-mer database class
    bp::class_<phylo_kmer_db, boost::noncopyable, std:: shared_ptr<phylo_kmer_db>>("PhyloKmerDb", bp::init<size_t, xpas::phylo_kmer::score_type, std::string>())
        .def("search", &phylo_kmer_db::search)
        .def("size", &phylo_kmer_db::size)
        .def("kmer_size", &phylo_kmer_db::kmer_size)
        .def("omega", &phylo_kmer_db::omega)

        .def("as_matrix", &db_as_matrix)

        .def("__iter__", bp::range(&phylo_kmer_db::begin, &phylo_kmer_db::end))
    ;

    /// database deserialization
    bp::def("load",
        +[](const std::string& filename)
        {
            /// core::phylo_kmer_db is move-able only. core::load returns an object, relying on
            /// guaranteed copy elision. But this mechanism does not exist in python.
            /// To avoid making phylo_kmer_db copyable, here we use smart pointers instead.
            return std::make_shared<phylo_kmer_db>(xpas::load(filename));
        });
}