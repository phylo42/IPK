#include <cmath>
#include <algorithm>
#include <xpas/seq.h>
#include <node_entry.h>
#include <window.h>
#include <cassert>

using namespace xpas;
using namespace xpas::impl;

kmer_iterator_pimpl::kmer_iterator_pimpl(window* window, size_t kmer_size, phylo_kmer::score_type threshold)
    : _window{window}
    , _kmer_size{kmer_size}
    , _threshold{threshold}
    , _current{ unpositioned_phylo_kmer::na_key, unpositioned_phylo_kmer::na_score }
{
}

kmer_iterator_pimpl& kmer_iterator_pimpl::operator=(kmer_iterator_pimpl&& rhs) noexcept
{
    if (*this != rhs)
    {
        _window = rhs._window;
        _kmer_size = rhs._kmer_size;
        _threshold = rhs._threshold;
        _current = rhs._current;
    }
    return *this;
}

bool kmer_iterator_pimpl::operator==(const kmer_iterator_pimpl& rhs) const noexcept
{
    return _window == rhs._window && _kmer_size == rhs._kmer_size && _threshold == rhs._threshold;
}

bool kmer_iterator_pimpl::operator!=(const kmer_iterator_pimpl& rhs) const noexcept
{
    return !(*this == rhs);
}

kmer_iterator_pimpl::reference kmer_iterator_pimpl::operator*() const noexcept
{
    return _current;
}

kmer_iterator_pimpl::pointer kmer_iterator_pimpl::operator->()  const noexcept
{
    return &_current;
}

///////////////////////////////////////////////////////////////////////////////

naive_DC_iterator xpas::impl::make_dac_end_iterator()
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
/// It does the same as dac_kmer_iterator::_next_kmer, but calculates all in one go, not on demand
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

naive_DC_iterator::naive_DC_iterator(window* window, size_t kmer_size, xpas::phylo_kmer::score_type threshold,
                                     xpas::phylo_kmer::pos_type start_pos, size_t prefix_size,
                                     vector_type<xpas::unpositioned_phylo_kmer> prefixes) noexcept
    : kmer_iterator_pimpl{ window, kmer_size, threshold }
    , _prefix_size{ 0 }
    , _start_pos{ start_pos }
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
        _prefixes = get_column(*_window->get_entry(), _start_pos);
        _prefix_it = _prefixes.begin();
        _current = _next_kmer();
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
                          : naive_DC_iterator(_window, _prefix_size, threshold, start_pos, 0, {});
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
                const auto central_column = get_column(*_window->get_entry(), _start_pos + _prefix_size);
                _prefixes = join_mmers(temp_prefixes, central_column, 1, _threshold);
                ++_prefix_size;
            }
            _prefix_it = _prefixes.begin();
        }

        // Calculate suffixes
        {
            auto it = (_prefix_size < kmer_size)
                      ? naive_DC_iterator(_window, kmer_size - _prefix_size, threshold, start_pos + _prefix_size, 0, {})
                      : make_dac_end_iterator();
            const auto end = make_dac_end_iterator();
            //std::cout << "\tsuffix of " << kmer_size - _prefix_size << " from " << start_pos + _prefix_size << ":" << std::endl;
            for (; it != end; ++it)
            {
                _suffixes.push_back(*it);
                //std::cout << "\t\t" << xpas::decode_kmer(it->key, kmer_size - _prefix_size)  << " " << std::pow(10, it->score) << std::endl;
            }

            _suffix_it = _suffixes.begin();
        }

        /// if there are no prefixes, the window is over
        if (_prefixes.empty())
        {
            _finish_iterator();
            _current = {};
        }
        else
        {
            _select_suffix_bound();
            _current = _next_kmer();
        }
    }
}

naive_DC_iterator& naive_DC_iterator::operator=(naive_DC_iterator&& rhs) noexcept
{
    if (*this != rhs)
    {
        _window = rhs._window;
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

naive_DC_iterator& naive_DC_iterator::operator++()
{
    _current = _next_kmer();
    return *this;
}

unpositioned_phylo_kmer naive_DC_iterator::_next_kmer()
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
                    << "= " << xpas::decode_kmer(full_key, _kmer_size)  << " " << std::pow(10, full_score) << std::endl;
                    */
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

void naive_DC_iterator::_select_suffix_bound()
{
    const auto residual_threshold = _threshold - _prefix_it->score;
    _last_suffix_it = std::partition(_suffixes.begin(), _suffixes.end(),
                                     [residual_threshold](auto pk) { return pk.score >= residual_threshold;});
}

void naive_DC_iterator::_finish_iterator()
{
    if (_window != nullptr)
    {
        /// Save suffixes of this window as prefixes of the next window
        _window->set_prefixes(std::move(_suffixes));
        /// save their length
        _window->set_prefix_size(_kmer_size - _prefix_size);
    }

    /// Now we can safely mark the iterator as ended
    *this = make_dac_end_iterator();
}

std::unique_ptr<kmer_iterator_pimpl> impl::make_iterator(const algorithm& algorithm, window* window, size_t start_pos,
                                                         size_t kmer_size, phylo_kmer::score_type threshold,
                                                         vector_type<xpas::unpositioned_phylo_kmer> prefixes)
{
    return std::make_unique<naive_DC_iterator>(window, kmer_size, threshold, start_pos, kmer_size / 2, std::move(prefixes));
}

kmer_iterator::kmer_iterator(const xpas::algorithm& algorithm, window* window, size_t kmer_size, phylo_kmer::score_type threshold)
    : _pimpl(impl::make_iterator(algorithm, window,
                                 (window == nullptr) ? 0 : (window->get_start_pos()), // for end() iterators there is no window
                                 kmer_size, threshold, {}))
{
}

kmer_iterator::reference kmer_iterator::operator*() const noexcept
{
    return _pimpl->operator*();
}

kmer_iterator::pointer kmer_iterator::operator->() const noexcept
{
    return _pimpl->operator->();
}

bool kmer_iterator::operator==(const kmer_iterator& rhs) const noexcept
{
    return *_pimpl == *(rhs._pimpl);
}

bool kmer_iterator::operator!=(const kmer_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

kmer_iterator& kmer_iterator::operator++()
{
    _pimpl->operator++();
    return *this;
}

enumerate_kmers::enumerate_kmers(const algorithm& algorithm, window* window,
                                 size_t kmer_size, phylo_kmer::score_type threshold,
                                 impl::vector_type<unpositioned_phylo_kmer> prefixes)
    : _algorithm(algorithm), _window(window), _kmer_size(kmer_size), _threshold(threshold)
{
}

enumerate_kmers::const_iterator enumerate_kmers::begin() const
{
    return { _algorithm, _window, _kmer_size, _threshold };
}

enumerate_kmers::const_iterator enumerate_kmers::end() const noexcept
{
    return { _algorithm, nullptr, 0, 0.0 };
}

window::window(const node_entry* entry, phylo_kmer::score_type threshold,
               phylo_kmer::pos_type start, phylo_kmer::pos_type end) noexcept
    : _entry{ entry }
    , _threshold{ threshold }
    , _start{ start }
    , _end{ end }
    /// _prefix_size will be set by the dac_iterator later
    , _prefix_size{ 0 }
    /// same for _prefixes. At the beginning it's empty
{}

window::window(const window& other) noexcept
    : window{ other._entry, other._threshold, other._start, other._end }
{
    /// We need the copy constructor, but only for windows without precomputed prefixes.
    /// Check that we never copy those around
    assert(other._prefixes.empty());
}


window& window::operator=(window&& other) noexcept
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

const node_entry* window::get_entry() const noexcept
{
    return _entry;
}

xpas::phylo_kmer::pos_type window::get_start_pos() const noexcept
{
    return _start;
}

void window::set_start_pos(xpas::phylo_kmer::pos_type pos)
{
    _start = pos;
}

xpas::phylo_kmer::pos_type window::get_end_pos() const noexcept
{
    return _end;
}

void window::set_end_pos(phylo_kmer::pos_type pos)
{
    _end = pos;
}

xpas::phylo_kmer::score_type window::get_threshold() const noexcept
{
    return _threshold;
}

void window::set_prefixes(impl::vector_type<xpas::unpositioned_phylo_kmer> prefixes)
{
    _prefixes = std::move(prefixes);
}

size_t window::get_prefix_size() const noexcept
{
    return _prefix_size;
}

void window::set_prefix_size(size_t size)
{
    _prefix_size = size;
}

bool xpas::operator==(const window& a, const window& b) noexcept
{
    return (a.get_start_pos() == b.get_start_pos()) &&
           (a.get_end_pos() == b.get_end_pos()) &&
           (a.get_entry() == b.get_entry());
}

bool xpas::operator!=(const window& a, const window& b) noexcept
{
    return !(a == b);
}
