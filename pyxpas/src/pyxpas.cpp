#include <boost/python.hpp>
#include <core/phylo_kmer_db.h>
#include <core/serialization.h>

using namespace core;
namespace bp = boost::python;

/// Custom exceptions
struct AttributeError: std::exception
{
    const char* what() const throw() { return "AttributeError exception"; }
};

struct TypeError: std::exception
{
    const char* what() const throw() { return "TypeError exception"; }
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

    static PyObject* convert(const std::optional<core::impl::search_result>& obj)
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
            core::impl::search_result empty_result(empty_collection.begin(), empty_collection.end());

            return bp::incref(bp::object(empty_result).ptr());
        }
    }
};

BOOST_PYTHON_MODULE(xpas)
{
    /// Exception translator
    bp::register_exception_translator<AttributeError>(&translate);
    bp::register_exception_translator<TypeError>(&translate);

    /// Python wrapping of std::optional.
    /// This is an adaption of the solution from stackoverlow.
    /// See https://stackoverflow.com/questions/36485840/wrap-boostoptional-using-boostpython
    bp::to_python_converter<std::optional<core::impl::search_result>, to_python_search_result>();

    /// a pair { branch -> score }
    bp::class_<pkdb_value>("PkDbValue", bp::init<phylo_kmer::branch_type, phylo_kmer::score_type>())
        .def_readonly("branch", &pkdb_value::branch)
        .def_readonly("score", &pkdb_value::score)
    ;

    /// a search result: collection of pkdb_value
    bp::class_<core::impl::search_result>("PkDbSearchResult")
        .add_property("values", bp::range(&core::impl::search_result::begin, &core::impl::search_result::end))
    ;

    /// phylo k-mer database class
    bp::class_<phylo_kmer_db, boost::noncopyable, std:: shared_ptr<phylo_kmer_db>>("PhyloKmerDb", bp::init<size_t, core::phylo_kmer::score_type, std::string>())
        .def("search", &phylo_kmer_db::search)
        .def("size", &phylo_kmer_db::size)
        .def("kmer_size", &phylo_kmer_db::kmer_size)
        .def("omega", &phylo_kmer_db::omega)
    ;

    /// database deserialization
    bp::def("load",
        +[](const std::string& filename)
        {
            /// core::phylo_kmer_db is move-able only. core::load returns an object, relying on
            /// guaranteed copy elision. But this mechanism does not exist in python.
            /// To avoid making phylo_kmer_db copyable, here we use smart pointers instead.
            return std::make_shared<phylo_kmer_db>(core::load(filename));
        });
}