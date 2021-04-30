#include <cmath>
#include <algorithm>
#include <xpas/seq.h>
#include <node_entry.h>
#include <node_entry_view.h>
#include <cassert>
#include <iostream>

using namespace xpas;
using namespace xpas::impl;


bnb_kmer_iterator make_bnb_end_iterator()
{
    return bnb_kmer_iterator();
}

bnb_kmer_iterator make_bnb_begin_iterator(const node_entry* entry, phylo_kmer::pos_type start, size_t kmer_size,
                                          phylo_kmer::score_type threshold)
{
    bnb_kmer_iterator::stack_type stack;
    stack.reserve(kmer_size + 1);

    /// Push a fake mmer to the bottom of the stack. We use it to reduce the number
    /// of if statements in the operator++, and therefore to reduce the number of branch
    /// mispredictions. It is also necessary to start branch-and-bound with every value
    /// of the first column without adding new if statements.
    stack.push_back(phylo_mmer{ { 0, 0.0 }, phylo_kmer::na_pos, 0, 0 });

    /// Calculate the first k-mer
    phylo_kmer::key_type kmer_key = 0;
    phylo_kmer::score_type kmer_score = 0.0;
    for (size_t i = 0; i < kmer_size; ++i)
    {
        const auto& ith_letter = entry->at(start + i, 0);
        kmer_key = (kmer_key << bit_length<seq_type>()) | ith_letter.index;
        kmer_score += ith_letter.score;

        stack.push_back(phylo_mmer{ { kmer_key, kmer_score }, phylo_kmer::pos_type(i), 0, 0 });
    }

    /// There is no point to iterate over the window if the first k-mer is already not good enough
    if (kmer_score < threshold)
    {
        return make_bnb_end_iterator();
    }
    else
    {
        return bnb_kmer_iterator{ entry, kmer_size, threshold, start, std::move(stack) };
    }
}

bnb_kmer_iterator::bnb_kmer_iterator() noexcept
    : _entry{ nullptr }, _kmer_size{ 0 }, _start_pos{ 0 }, _threshold{ 0.0 }, _stack{ }
{}

bnb_kmer_iterator::bnb_kmer_iterator(const node_entry* entry, size_t kmer_size, phylo_kmer::score_type threshold,
                                     phylo_kmer::pos_type start_pos, stack_type&& stack) noexcept
    : _entry{ entry }
      , _kmer_size{ kmer_size }
      , _start_pos{ start_pos }
      , _threshold{ threshold }
      , _stack{ std::move(stack) }
{
    /// stack can be empty for the end() method
    if (!_stack.empty())
    {
        _current = _stack.back();
    }
}

bnb_kmer_iterator& bnb_kmer_iterator::operator=(bnb_kmer_iterator&& rhs) noexcept
{
    if (*this != rhs)
    {
        _entry = rhs._entry;
        _kmer_size = rhs._kmer_size;
        _start_pos = rhs._start_pos;
        _threshold = rhs._threshold;
        _stack = std::move(rhs._stack);
        _current = rhs._current;
    }
    return *this;
}

bool bnb_kmer_iterator::operator==(const bnb_kmer_iterator& rhs) const noexcept
{
    if (_entry != rhs._entry)
    {
        return false;
    }
    if (_stack.empty())
    {
        return rhs._stack.empty();
    }
    if (rhs._stack.empty())
    {
        return false;
    }
    return _stack.back().mmer == rhs._stack.back().mmer;
}

bool bnb_kmer_iterator::operator!=(const bnb_kmer_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

bnb_kmer_iterator& bnb_kmer_iterator::operator++()
{
    _current = next_phylokmer();
    return *this;
}

bnb_kmer_iterator::reference bnb_kmer_iterator::operator*() const noexcept
{
    return _current.mmer;
}

bnb_kmer_iterator::pointer bnb_kmer_iterator::operator->() const noexcept
{
    return &(_current.mmer);
}

phylo_mmer bnb_kmer_iterator::next_phylokmer()
{
    {
        const auto next_index = _stack.back().last_index + 1;
        _stack.pop_back();
        _stack.back().next_index = next_index;
    }

    /// until there is only the fake k-mer
    while (!(_stack.size() == 1 && _stack.back().next_index == seq_traits::alphabet_size))
    {
        const auto top_mmer = _stack.back();

        /// if we found a m-mer with low score, we cut the whole subtree
        if (top_mmer.mmer.score < _threshold)
        {
            _stack.pop_back();
            _stack.back().next_index = seq_traits::alphabet_size;
        }
        else
        {
            /// |top_mmer| < k, go to the next column
            /// + 1 because of the fake mmer
            if (_stack.size() < _kmer_size + 1)
            {
                /// go to the next letter in the next column
                if (top_mmer.next_index < seq_traits::alphabet_size)
                {
                    const phylo_kmer::pos_type new_letter_position = top_mmer.last_position + 1;
                    const auto& new_letter = _entry->at(_start_pos + new_letter_position, top_mmer.next_index);
                    const auto new_mmer_key = (top_mmer.mmer.key << bit_length<seq_type>()) | new_letter.index;
                    const auto new_mmer_score = top_mmer.mmer.score + new_letter.score;
                    _stack.push_back(phylo_mmer{ { new_mmer_key, new_mmer_score },
                                                 new_letter_position, top_mmer.next_index, 0 });
                }
                /// the next column is over, get back
                else
                {
                    _stack.pop_back();
                    _stack.back().next_index = top_mmer.last_index + 1;
                }
            }
            /// we have a good k-mer
            else
            {
                return top_mmer;
            }
        }
    }

    *this = make_bnb_end_iterator();
    return {};
}


dac_kmer_iterator xpas::impl::make_dac_end_iterator()
{
    return { nullptr, 0, 0.0, 0, 0, {} };
}

bool kmer_score_comparator(const unpositioned_phylo_kmer& k1, const unpositioned_phylo_kmer& k2)
{
    return k1.score > k2.score;
}

bool is_odd(int x)
{
    return x % 2;
}

/// Creates a vector of 1-mers from a column of PP matrix
vector_type<xpas::unpositioned_phylo_kmer> get_column(const node_entry& entry, size_t pos)
{
    vector_type<xpas::unpositioned_phylo_kmer> column;
    for (size_t i = 0; i < seq_traits::alphabet_size; ++i)
    {
        const auto& letter = entry.at(pos, i);
        column.push_back(make_phylo_kmer<unpositioned_phylo_kmer>(letter.index, letter.score, pos));
    }
    return column;
}

/// Creates all combinations of prefixes and suffixes
/// It does the same as dac_kmer_iterator::_next_phylokmer, but calculates all in one go, not on demand
vector_type<xpas::unpositioned_phylo_kmer> join_mmers(const vector_type<xpas::unpositioned_phylo_kmer>& prefixes,
                                                      const vector_type<xpas::unpositioned_phylo_kmer>& suffixes,
                                                      size_t suffix_size,
                                                      phylo_kmer::score_type threshold)
{
    vector_type<xpas::unpositioned_phylo_kmer> result;

    for (const auto& prefix : prefixes)
    {
        const auto residual_threshold = threshold - prefix.score;

        for (const auto& suffix : suffixes)
        {
            if (suffix.score < residual_threshold)
            {
                break;
            }

            const auto full_key = (prefix.key << (suffix_size * xpas::bit_length<xpas::seq_type>())) | suffix.key;
            const auto full_score = prefix.score + suffix.score;
            result.push_back(make_phylo_kmer<unpositioned_phylo_kmer>(full_key, full_score, 0));
        }
    }
    return result;
}

dac_kmer_iterator::dac_kmer_iterator(node_entry_view* view, size_t kmer_size, xpas::phylo_kmer::score_type threshold,
                                     xpas::phylo_kmer::pos_type start_pos, size_t prefix_size,
                                     vector_type<xpas::unpositioned_phylo_kmer> prefixes) noexcept
    : _entry_view{ view }, _kmer_size{ kmer_size }, _prefix_size{ 0 }, _start_pos{ start_pos }, _threshold{ threshold }
{
    /// if prefixes are not given, figure out their size
    if (prefixes.empty())
    {
        /// kmer_size can also be zero, which means the end() iterator
        const auto halfsize = size_t{ kmer_size / 2 };
        _prefix_size = (halfsize >= 1) ? halfsize : kmer_size;
    }
    else
    {
        _prefix_size = prefix_size;
    }

    /// for the window of size 1, k-mers are trivial
    if (kmer_size == 1)
    {
        _prefixes = get_column(*_entry_view->get_entry(), _start_pos);
        _prefix_it = _prefixes.begin();
        _current = _next_phylokmer();
    }
    /// for all other window sizes:
    /// if it is not the end()
    else if (_prefix_size > 0)
    {
        // Calculate prefixes
        {
            /// if prefixes are provided, they are already calculated from the previous window
            if (!prefixes.empty())
            {
                _prefixes = std::move(prefixes);
            }
            /// otherwise it's the very first window in the chain: calculate prefixes
            else
            {
                auto it = (_prefix_size == 0)
                          ? make_dac_end_iterator()
                          : dac_kmer_iterator(_entry_view, _prefix_size, threshold, start_pos, 0, {});
                const auto end = make_dac_end_iterator();
                //std::cout << "\tprefix of " << _prefix_size << " from " << start_pos << ":" << std::endl;
                for (; it != end; ++it)
                {
                    _prefixes.push_back(*it);
                    //std::cout << "\t\t" << xpas::decode_kmer(it->key, _prefix_size)  << " " << std::pow(10, it->score) << std::endl;
                }
            }

            /// if k is odd, we need to join one more column to the window of prefixes
            if (is_odd(kmer_size))
            {
                const auto temp_prefixes = std::move(_prefixes);
                const auto central_column = get_column(*_entry_view->get_entry(), _start_pos + _prefix_size);
                _prefixes = join_mmers(temp_prefixes, central_column, 1, _threshold);
                ++_prefix_size;
            }
            _prefix_it = _prefixes.begin();
        }

        /// if there are no prefixes, the window is over
        if (_prefixes.empty())
        {
            _suffixes.clear();
            _finish_iterator();
            _current = {};
        }
        // Calculate suffixes
        else
        {
            auto it = (_prefix_size < kmer_size)
                      ? dac_kmer_iterator(_entry_view, kmer_size - _prefix_size, threshold, start_pos + _prefix_size, 0, {})
                      : make_dac_end_iterator();
            const auto end = make_dac_end_iterator();
            //std::cout << "\tsuffix of " << kmer_size - _prefix_size << " from " << start_pos + _prefix_size << ":" << std::endl;
            for (; it != end; ++it)
            {
                _suffixes.push_back(*it);
                //std::cout << "\t\t" << xpas::decode_kmer(it->key, kmer_size - _prefix_size)  << " " << std::pow(10, it->score) << std::endl;
            }

            _suffix_it = _suffixes.begin();
            _select_suffix_bound();
            _current = _next_phylokmer();
        }
    }
}

dac_kmer_iterator& dac_kmer_iterator::operator=(dac_kmer_iterator&& rhs) noexcept
{
    if (*this != rhs)
    {
        _entry_view = rhs._entry_view;
        _kmer_size = rhs._kmer_size;
        _prefix_size = rhs._prefix_size;
        _start_pos = rhs._start_pos;
        _threshold = rhs._threshold;
        _current = rhs._current;
        _prefixes = std::move(rhs._prefixes);
        _prefix_it = rhs._prefix_it;
        _suffixes = std::move(rhs._suffixes);
        _suffix_it = rhs._suffix_it;
        _last_suffix_it = rhs._last_suffix_it;
    }
    return *this;
}

bool dac_kmer_iterator::operator==(const dac_kmer_iterator& rhs) const noexcept
{
    return _entry_view == rhs._entry_view && _start_pos == rhs._start_pos && _kmer_size == rhs._kmer_size;
}

bool dac_kmer_iterator::operator!=(const dac_kmer_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

dac_kmer_iterator& dac_kmer_iterator::operator++()
{
    _current = _next_phylokmer();
    return *this;
}

dac_kmer_iterator::reference dac_kmer_iterator::operator*() const noexcept
{
    return _current;
}

dac_kmer_iterator::pointer dac_kmer_iterator::operator->()  const noexcept
{
    return &_current;
}

unpositioned_phylo_kmer dac_kmer_iterator::_next_phylokmer()
{
    if (_prefix_it != _prefixes.end())
    {
        if (_kmer_size > 1)
        {
            while ((_suffix_it == _last_suffix_it) && (_prefix_it != _prefixes.end()))
            {
                ++_prefix_it;
                _suffix_it = _suffixes.begin();
                _select_suffix_bound();
            }

            if (_prefix_it != _prefixes.end())
            {
                const auto prefix = *_prefix_it;
                const auto suffix = *_suffix_it;
                const auto full_key =
                    (prefix.key << ((_kmer_size - _prefix_size) * xpas::bit_length<xpas::seq_type>()))
                    | suffix.key;
                const auto full_score = prefix.score + suffix.score;
                ++_suffix_it;
                /*std::cout << "\t\tpairing "
                    << xpas::decode_kmer(prefix.key, _prefix_size)  << " " << std::pow(10, prefix.score) << " "
                    << xpas::decode_kmer(suffix.key, _kmer_size - _prefix_size)  << " " << std::pow(10, suffix.score) << " "
                    << "= " << xpas::decode_kmer(full_key, _kmer_size)  << " " << std::pow(10, full_score) << std::endl;*/
                return make_phylo_kmer<unpositioned_phylo_kmer>(full_key, full_score, 0);
            }
        }
        else
        {
            const auto kmer = *_prefix_it;
            ++_prefix_it;
            return kmer;
        }
    }

    /// If we reached this line, the iterator is over
    _finish_iterator();
    return {};
}

void dac_kmer_iterator::_select_suffix_bound()
{
    const auto residual_threshold = _threshold - _prefix_it->score;
    _last_suffix_it = ::std::lower_bound(_suffixes.begin(), _suffixes.end(),
                                         make_phylo_kmer<unpositioned_phylo_kmer>(0, residual_threshold, 0), kmer_score_comparator);
}

void dac_kmer_iterator::_finish_iterator()
{
    if (_entry_view != nullptr)
    {
        /// Save suffixes of this window as prefixes of the next window
        _entry_view->set_prefixes(std::move(_suffixes));
        /// save their length
        _entry_view->set_prefix_size(_kmer_size - _prefix_size);
    }

    /// Now we can safely mark the iterator as ended
    *this = make_dac_end_iterator();
}

node_entry_view::node_entry_view(const node_entry* entry, phylo_kmer::score_type threshold,
                                 phylo_kmer::pos_type start, phylo_kmer::pos_type end) noexcept
    : _entry{ entry }
    , _threshold{ threshold }
    , _start{ start }
    , _end{ end }
    /// _prefix_size will be set by the dac_iterator later
    , _prefix_size{ 0 }
    /// same for _prefixes. At the beginning it's empty
{}

node_entry_view::node_entry_view(const node_entry_view& other) noexcept
    : node_entry_view{ other._entry, other._threshold, other._start, other._end }
{
    /// We need the copy constructor, but only for windows without precomputed prefixes.
    /// Check that we never copy those around
    assert(other._prefixes.empty());
}

node_entry_view::iterator node_entry_view::begin()
{
    const auto kmer_size = size_t{ (size_t)_end - _start + 1};
    //const node_entry* entry, size_t kmer_size, phylo_kmer::score_type threshold,
    //                                     phylo_kmer::pos_type start_pos, stack_type&& stack)

    // DAC-CW:
    //return { this, kmer_size, _threshold, _start, _prefix_size, std::move(_prefixes) };

    // DAC:
    return { this, kmer_size, _threshold, _start, _prefix_size, {} };

    // BNB:
    //return make_bnb_begin_iterator(_entry, _start, kmer_size, _threshold);
}

node_entry_view::iterator node_entry_view::end() const noexcept
{
    // DAC:
    return make_dac_end_iterator();
    // BNB:
    //return make_bnb_end_iterator();
}

node_entry_view& node_entry_view::operator=(node_entry_view&& other) noexcept
{
    if (*this != other)
    {
        _entry = other._entry;
        _start = other._start;
        _end = other._end;
        other._entry = nullptr;
        other._start = 0;
        other._end = 0;
    }
    return *this;
}

const node_entry* node_entry_view::get_entry() const noexcept
{
    return _entry;
}

xpas::phylo_kmer::pos_type node_entry_view::get_start_pos() const noexcept
{
    return _start;
}

void node_entry_view::set_start_pos(xpas::phylo_kmer::pos_type pos)
{
    _start = pos;
}

xpas::phylo_kmer::pos_type node_entry_view::get_end_pos() const noexcept
{
    return _end;
}

void node_entry_view::set_end_pos(phylo_kmer::pos_type pos)
{
    _end = pos;
}

xpas::phylo_kmer::score_type node_entry_view::get_threshold() const noexcept
{
    return _threshold;
}

void node_entry_view::set_prefixes(impl::vector_type<xpas::unpositioned_phylo_kmer> prefixes)
{
    _prefixes = std::move(prefixes);
}

size_t node_entry_view::get_prefix_size() const noexcept
{
    return _prefix_size;
}

void node_entry_view::set_prefix_size(size_t size)
{
    _prefix_size = size;
}

bool xpas::operator==(const node_entry_view& a, const node_entry_view& b) noexcept
{
    return (a.get_start_pos() == b.get_start_pos()) &&
           (a.get_end_pos() == b.get_end_pos()) &&
           (a.get_entry() == b.get_entry());
}

bool xpas::operator!=(const node_entry_view& a, const node_entry_view& b) noexcept
{
    return !(a == b);
}
