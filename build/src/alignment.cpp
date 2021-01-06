#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <boost/filesystem.hpp>

#include "fasta.h"
#include "alignment.h"
#include "phylo_tree.h"

using std::vector;
using std::string;
using namespace xpas;
namespace fs = boost::filesystem;

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


//------------------------------------------------------------------------------------
// alignment class
alignment::alignment(std::vector<seq_record> sequences)
    : _sequences(std::move(sequences))
{
    if (_sequences.empty())
    {
        throw std::runtime_error("The alignment is empty.");
    }
    _width = _sequences[0].sequence().size();
    _height = _sequences.size();
}

size_t alignment::width() const noexcept
{
    return _width;
}

size_t alignment::height() const noexcept
{
    return _sequences.size();
}

alignment::iterator alignment::begin()
{
    return _sequences.begin();
}

alignment::iterator alignment::end()
{
    return _sequences.end();
}

alignment::const_iterator alignment::begin() const
{
    return _sequences.begin();
}

alignment::const_iterator alignment::end() const
{
    return _sequences.end();
}

//------------------------------------------------------------------------------------
// miscellaneous functions
alignment load_alignment(const string& file_name)
{
    std::vector<seq_record> sequences;
    for (const auto& seq : xpas::io::read_fasta(file_name, false))
    {
        sequences.push_back(seq);
    }

    return alignment(std::move(sequences));
}

// save a stream of fasta records to a file
template <class InputIt>
void save_fasta(InputIt first, InputIt last, const std::string& file_name)
{

    std::ofstream out(file_name);
    for (auto it = first; it != last; ++it)
    {
        out << ">" << it->header() << '\n' << it->sequence() << '\n';
    }
}

// save a stream of fasta records to a file
template <class InputIt>
void save_phylip(InputIt first, InputIt last, const std::string& file_name)
{
    const size_t allowed_label_size = 250;

    std::ofstream out(file_name);

    /// Header
    out << '\t' << std::distance(first, last) << '\t' << first->sequence().size() << '\n';

    for (auto it = first; it != last; ++it)
    {
        /// Left column
        out << it->header();
        for (size_t i = it->header().size(); i < allowed_label_size; ++i)
        {
            out << ' ';
        }

        /// Other columns
        size_t pos = 0;
        while (pos < it->sequence().size())
        {
            size_t remained = it->sequence().size() - pos;
            if (remained > 10)
            {
                out << it->sequence().substr(pos, 10) << ' ';
                pos += 10;
            }
            else
            {
                out << it->sequence().substr(pos, remained);
                pos += remained;
            }
        }
        out << '\n';

    }
}

void xpas::save_alignment(const alignment& align, const string& file_name, alignment_format format)
{
    if (format == alignment_format::FASTA)
    {
        save_fasta(std::begin(align), std::end(align), file_name);
    }
    else if (format == alignment_format::PHYLIP)
    {
        save_phylip(std::begin(align), std::end(align), file_name);
    }
}

vector<double> calculate_gap_ratio(const alignment& align)
{
    auto ratios = vector<double>(align.width(), 0.0f);

    for (const auto& seq_record : align)
    {
        const auto& sequence = seq_record.sequence();
        for (size_t i = 0; i < sequence.size(); ++i)
        {
            if (xpas::seq_traits::is_gap(sequence[i]))
            {
                ratios[i]++;
            }
        }
    }

    for (double& ratio : ratios)
    {
        ratio /= (double)align.height();
    }
    return ratios;
}

alignment reduce_alignment(const alignment &align, double reduction_ratio)
{
    // figure out which columns should be removed
    const auto gap_ratios = calculate_gap_ratio(align);
    vector<bool> pos_to_remove(gap_ratios.size(), false);
    transform(begin(gap_ratios), end(gap_ratios), pos_to_remove.begin(),
              [reduction_ratio](double ratio) -> bool { return ratio >= reduction_ratio; });

    vector<seq_record> reduced_sequences;
    for (const auto& seq_record : align)
    {
        auto header = std::string(seq_record.header());
        auto sequence = std::string(seq_record.sequence());

        // erase-remove idiom: remove the positions, marked as 1 in the pos_to_remove vector
        size_t index = 0;
        auto predicate = [&index, pos_to_remove](char) -> bool { return pos_to_remove[index++]; };
        sequence.erase(
            remove_if(sequence.begin(), sequence.end(), predicate),
            sequence.end()
        );

        reduced_sequences.emplace_back(move(header), move(sequence));
    }
    return alignment(move(reduced_sequences));
}

void check_length(const alignment& align)
{
    for (const auto &sequence : align)
    {
        if (sequence.sequence().size() != align.width())
        {
            const auto header = std::string(sequence.header());
            const auto first_header = std::string(std::begin(align)->header());
            std::ostringstream ss;
            ss << "Error: Sequences in the input alignment do not have same number of sites. "
               << header << " is " << sequence.sequence().size() << "bp in length, while "
               << first_header << " is " << align.width() << "bp in length.";
            throw std::runtime_error(ss.str());
        }
    }
}

void check_sequence_states(const seq_record& seq_record)
{
    for (const auto& state : seq_record.sequence())
    {
        if (!xpas::seq_traits::is_gap(state) && !xpas::seq_traits::key_to_code(state))
        {
            std::ostringstream ss;
            ss << "Error: " << std::string(seq_record.header()) << " contains a non supported state: " <<
                state;
            throw std::runtime_error(ss.str());
        }
        else if (xpas::seq_traits::is_ambiguous(state))
        {
            if (!xpas::seq_traits::is_gap(state))
            {
                // TODO: logging
                //std::cout << "Ambiguous state (char='" << state << "') will be considered as a gap during AR." << std::endl;
            }
        }
    }
}

void check_sequence_states(const alignment& align)
{
    for (const auto& sequence : align)
    {
        check_sequence_states(sequence);
    }
}

void validate_alignment(const alignment& align)
{
    // test if the sequence lengths are the same
    check_length(align);

    // test if the sequences have unsupported states
    //check_sequence_states(align);
}

xpas::alignment _preprocess_alignment(const std::string& working_dir,
                                      const std::string& alignment_file,
                                      double reduction_ratio,
                                      bool no_reduction)
{
    /// Create working directories
    fs::create_directories(working_dir);

    /// Read and process the reference alignment
    auto raw_alignment = load_alignment(alignment_file);
    validate_alignment(raw_alignment);

    /// Reduce the alignment if needed
    if (no_reduction)
    {
        return raw_alignment;
    }
    else
    {
        auto alignment = reduce_alignment(raw_alignment, reduction_ratio);
        validate_alignment(alignment);

        /// Save the reduced alignment on disk
        const auto reduced_alignment_file = (fs::path(working_dir) / "align.reduced.fasta").string();
        save_alignment(alignment, reduced_alignment_file, alignment_format::FASTA);

        return alignment;
    }

}

xpas::alignment xpas::preprocess_alignment(const std::string& working_dir,
                                           const std::string& alignment_file,
                                           double reduction_ratio,
                                           bool no_reduction)
{
    std::cout << "Loading the reference alignment: " << alignment_file <<  std::endl;
    auto alignment = _preprocess_alignment(working_dir, alignment_file, reduction_ratio, no_reduction);
    std::cout << "Loaded and filtered " << alignment.height() <<  " sequences." << std::endl << std::endl;
    return alignment;
}



bool has_sequence(const alignment& alignment, const std::string& seq_header)
{
    return std::find_if(alignment.begin(), alignment.end(), [&seq_header](const seq_record& seq_rec) {
        return seq_rec.header() == seq_header;
    }) != alignment.end();
}

alignment xpas::extend_alignment(alignment original_alignment, const phylo_tree& tree)
{
    alignment extended_alignment = std::move(original_alignment);

    std::string empty_seq(extended_alignment.width(), seq_traits::get_gap());

    for (const auto& node : visit_subtree(tree.get_root()))
    {
        if (node.is_leaf() && !has_sequence(extended_alignment, node.get_label()))
        {
            extended_alignment._sequences.emplace_back(node.get_label(), empty_seq);
            extended_alignment._height += 1;
        }
    }

    return extended_alignment;
}