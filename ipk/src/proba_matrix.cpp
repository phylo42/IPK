#include "proba_matrix.h"
#include "ar.h"

using namespace ipk;

proba_matrix::proba_matrix(std::unique_ptr<ipk::ar::reader> reader)
    : _reader(std::move(reader))
{
}

size_t proba_matrix::num_branches() const
{
    return _data.size();
}

size_t proba_matrix::num_sites() const
{
    return std::begin(_data)->second.width();
}

proba_matrix::mapped_type& proba_matrix::operator[](const std::string& ar_label)
{
    return _data[ar_label];
}

const proba_matrix::mapped_type& proba_matrix::at(const std::string& ar_label) const
{
    return _data.at(ar_label);
}

proba_matrix::iterator proba_matrix::find(const std::string& ar_label)
{
    auto it = _data.find(ar_label);
    if (it == _data.end())
    {
        _data[ar_label] = _reader->read_node(ar_label);
        it = _data.find(ar_label);
    }
    return it;
}

proba_matrix::iterator proba_matrix::begin()
{
    return std::begin(_data);
}

proba_matrix::iterator proba_matrix::end()
{
    return std::end(_data);
}

proba_matrix::const_iterator proba_matrix::begin() const
{
    return _data.cbegin();
}

proba_matrix::const_iterator proba_matrix::end() const
{
    return std::end(_data);
}
