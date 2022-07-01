#include <xcl/seq_record.h>

using namespace xpas;
using std::string;

//------------------------------------------------------------------------------------
// seq_record
seq_record::seq_record(string header, string sequence) noexcept
    : _header(move(header)), _sequence(move(sequence))
{
}

std::string_view seq_record::header() const noexcept
{
    return _header;
}

std::string_view seq_record::sequence() const noexcept
{
    return _sequence;
}

bool seq_record::operator==(const seq_record& rhs) const noexcept
{
    return _header == rhs.header() && _sequence == rhs.sequence();
}

