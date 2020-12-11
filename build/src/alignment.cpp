#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <xpas/fasta.h>

#include "alignment.h"
#include "io.h"

using std::vector;
using std::string;
using xpas::io::fasta;
using namespace rappas;
namespace fs = boost::filesystem;

//------------------------------------------------------------------------------------
// alignment class
alignment::alignment(std::vector<fasta> sequences)
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
    std::vector<fasta> sequences;
    for (const auto& seq : xpas::io::read_fasta(file_name))
    {
        sequences.push_back(seq);
    }

    return alignment(std::move(sequences));
}

// save a stream of fasta records to a file
template <class InputIt>
void save_fasta(InputIt first, InputIt last, const std::string& file_name)
{
    std::cout << "Saving alignment to " << file_name << "..." << std::endl;

    std::ofstream out(file_name);
    for (auto it = first; it != last; ++it)
    {
        out << ">" << it->header() << std::endl << it->sequence() << std::endl;
    }
}

void save_alignment(const alignment& align, const string& file_name)
{
    save_fasta(std::begin(align), std::end(align), file_name);
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

    vector<fasta> reduced_sequences;
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

void check_sequence_states(const fasta& fasta)
{
    for (const auto& state : fasta.sequence())
    {
        if (!xpas::seq_traits::is_gap(state) && !xpas::seq_traits::key_to_code(state))
        {
            std::ostringstream ss;
            ss << "Error: " << std::string(fasta.header()) << " contains a non supported state: " <<
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

rappas::alignment _preprocess_alignment(const std::string& working_dir,
                                        const std::string& alignment_file,
                                        double reduction_ratio,
                                        bool no_reduction)
{
    /// Create working directories
    rappas::io::create_directory(working_dir);

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
        const auto reduced_alignment_file = (fs::path(working_dir) / "align.reduced").string();
        save_alignment(alignment, reduced_alignment_file);

        return alignment;
    }

}

rappas::alignment rappas::preprocess_alignment(const std::string& working_dir,
                                               const std::string& alignment_file,
                                               double reduction_ratio,
                                               bool no_reduction)
{
    std::cout << "Loading the reference alignment: " << alignment_file <<  std::endl;
    auto alignment = _preprocess_alignment(working_dir, alignment_file, reduction_ratio, no_reduction);
    std::cout << "Loaded and filtered " << alignment.height() <<  " sequences." << std::endl << std::endl;
    return alignment;
}