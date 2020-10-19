#include <iostream>
#include <boost/algorithm/string/predicate.hpp>
#include <utils/io/fasta.h>

using namespace xpas::io;
using std::vector;
using std::string, std::string_view;
using std::cout, std::endl;
using std::move;


fasta::fasta(string&& header, string&& sequence) noexcept
    : _header(move(header)), _sequence(move(sequence))
{
}

string_view fasta::header() const noexcept
{
    return _header;
}

string_view fasta::sequence() const noexcept
{
    return _sequence;
}

bool fasta::operator==(const fasta& rhs) const noexcept
{
    return _header == rhs.header() && _sequence == rhs.sequence();
}

//------------------------------------------------------------------------------------

impl::fasta_iterator::fasta_iterator(const std::string& filename, size_t batch_size)
    : _mmap{ filename }
    , _is{ _mmap, std::ios::in }
    , _batch_size{ batch_size }
    , _seq_id{ 0 }
    , _global_seq_id{ 0 }
    , _last_batch{ false }
{
    /// end() iterator
    if (_batch_size == 0)
    {
        _global_seq_id = 0;
    }
    else
    {
        _read_batch();
    }
}

impl::fasta_iterator& impl::fasta_iterator::operator++()
{
    ++_seq_id;
    ++_global_seq_id;

    /// if end of file, make it an end iterator
    if (_last_batch && _seq_id == _batch_size)
    {
        _global_seq_id = 0;
        _batch_size = 0;
    }
    else
    {
        if (_seq_id == _batch_size)
        {
            _read_batch();
        }
    }
    return *this;
}

bool impl::fasta_iterator::operator==(const fasta_iterator& rhs) const noexcept
{
    return _global_seq_id == rhs._global_seq_id && _batch_size == rhs._batch_size;
}

bool impl::fasta_iterator::operator!=(const fasta_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

impl::fasta_iterator::value_type impl::fasta_iterator::operator*() const noexcept
{
    return _seqs[_seq_id];
}

void impl::fasta_iterator::_read_batch()
{
    /// start iterating from the first sequence of the new batch
    _seq_id = 0;

    /// clean the old batch
    _seqs.clear();
    _seqs.reserve(_batch_size);

    /// read sequences until eof or batch limit
    string line, sequence;
    size_t i = 0;
    while ((i < _batch_size) && std::getline(_is, line))
    {
        if (boost::starts_with(line, ">"))
        {
            if (!sequence.empty())
            {
                ++i;
                _seqs.emplace_back(move(_header), clean_sequence(sequence));
                sequence = "";
                sequence.reserve(1024);
            }

            _header = {line.c_str() + 1, line.size() - 1};
        }
        else if (!line.empty())
        {
            sequence.append(line);
        }
    }

    /// if eof, it was the last batch
    if (i < _batch_size)
    {
        /// do not forget the last sequence
        _seqs.emplace_back(move(_header), clean_sequence(sequence));
        _last_batch = true;
        _batch_size = i + 1;
    }
}

read_fasta::read_fasta(std::string filename, size_t batch_size)
    : _filename{ std::move(filename) }
    , _batch_size{ batch_size }
{}

read_fasta::const_iterator read_fasta::begin() const
{
    return { _filename, _batch_size };
}

read_fasta::const_iterator read_fasta::end() const
{
    return { _filename, 0 };
}

string xpas::io::clean_sequence(string sequence)
{
    sequence.erase(
        std::remove_if(sequence.begin(), sequence.end(),
                       [](unsigned char x){ return x == '-' || x == '.' || x == '!' || x == '*'; }),
        sequence.end());
    return sequence;
}


