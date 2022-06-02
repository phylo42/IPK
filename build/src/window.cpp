#include <cmath>
#include <algorithm>
#include <xpas/seq.h>
#include <window.h>
#include <cassert>

using namespace xpas;
using namespace xpas::impl;

matrix::matrix(std::vector<column> data, std::string label)
    : _data(std::move(data)), _label(std::move(label))
{
    preprocess();
}

void matrix::preprocess()
{
    _best_scores = std::vector<score_t>(_data.size() + 1, 0.0f);
    score_t product = 0.0f;
    for (size_t j = 0; j < _data.size(); ++j)
    {
        const auto& column = get_column(j);
        const auto best_score = *std::max_element(column.begin(), column.end());
        product += best_score;
        _best_scores[j + 1] = product;
    }
}

score_t matrix::get(size_t i, size_t j) const
{
    return _data[j][i];
}

size_t matrix::width() const
{
    return _data.size();
}

bool matrix::empty() const
{
    return _data.empty();
}

void matrix::set_label(const std::string& label)
{
    _label = label;
}

const std::vector<matrix::column>& matrix::get_data() const
{
    return _data;
}

std::vector<matrix::column>& matrix::get_data()
{
    return _data;
}

std::string matrix::get_label() const
{
    return _label;
}

const matrix::column& matrix::get_column(size_t j) const
{
    return _data[j];
}

score_t matrix::range_max_sum(size_t start_pos, size_t len) const
{
    return _best_scores[start_pos + len] - _best_scores[start_pos];
}

window::window(const matrix* m, size_t start_pos, size_t size)
    : _matrix(m), _start_pos(start_pos), _size(size)
{
}
/*
window& window::operator=(window&& other) noexcept
{
    if (*this != other)
    {
        _matrix = std::move(other._matrix);
        _start_pos = other._start_pos;
        _size = other._size;
        _best_scores = other._best_scores;
    }
    return *this;
}*/

bool window::operator==(const window& other) const
{
    /// FIXME:
    /// Let us imagine those windows are over the same matrix
    return _start_pos == other._start_pos && _size == other._size;
}

bool window::operator!=(const window& other) const
{
    return !(*this == other);
}

score_t window::get(size_t i, size_t j) const
{
    return _matrix->get(i, _start_pos + j);
}

size_t window::size() const
{
    return _size;
}

bool window::empty() const
{
    return _size == 0;
}

size_t window::get_position() const
{
    return _start_pos;
}

phylo_kmer::score_type window::range_max_product(size_t pos, size_t len) const
{
    return _matrix->range_max_sum(_start_pos + pos, len);
}

matrix::column window::get_column(size_t j) const
{
    return _matrix->get_column(j);
}

std::pair<size_t, score_t> window::max_at(size_t column) const
{
    size_t max_index = 0;
    score_t max_score = get(0, column);
    for (size_t i = 1; i < seq_traits::alphabet_size; ++i)
    {
        if (get(i, column) > max_score)
        {
            max_score = get(i, column);
            max_index = i;
        }
    }
    return { max_index, max_score };
}

impl::window_iterator::window_iterator(const matrix* matrix, size_t kmer_size) noexcept
    : _matrix(matrix), _window(matrix, 0, kmer_size), _kmer_size(kmer_size), _current_pos(0)
{
}

impl::window_iterator& impl::window_iterator::operator++()
{
    _current_pos++;
    if (_current_pos + _kmer_size <= _matrix->width())
    {
        //_window = window(_matrix, _current_pos, _kmer_size);
        _window._start_pos = _current_pos;
        _window._size = _kmer_size;
    }
    else
    {
        // end iterator
        //_window = window(_matrix, 0, 0);
        _window._start_pos = 0;
        _window._size = 0;
        _kmer_size = 0;
    }
    return *this;
}

bool impl::window_iterator::operator==(const window_iterator& rhs) const noexcept
{
    return _window == rhs._window;
}

bool impl::window_iterator::operator!=(const window_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

impl::window_iterator::reference impl::window_iterator::operator*() noexcept
{
    return _window;
}


impl::chained_window_iterator::chained_window_iterator(const matrix* matrix, size_t kmer_size)
    : _matrix(matrix),
    _window(matrix, 0, kmer_size),
    _previous_window(matrix, 0, 0),
    _next_window(matrix, 0, 0),
    _kmer_size(kmer_size)
{
    if (_kmer_size > matrix->width())
    {
        throw std::runtime_error("Window is too small");
    }

    /// The start position of the last possible chain
    /// FIXME: odd k only
    _last_chain_pos = _kmer_size / 2 - 1;
    _chain_start = 0;

    // There is no previous window
    _previous_window = window(_matrix, 0, 0);

    const auto prefix_size = _kmer_size / 2;
    const auto suffix_size = _kmer_size - prefix_size;
    /// see if the next window is in the current chain
    if (size_t(_window.get_position() + suffix_size + _kmer_size) <= _matrix->width())
    {
        _next_window = { _matrix, _window.get_position() + suffix_size, _kmer_size };
    }
    /// if not, see if there is another chain
    else if ((_chain_start + 1 <= _last_chain_pos) && (_chain_start + 1 + _kmer_size <= _matrix->width()))
    {
        _next_window = { _matrix, _chain_start + 1, _kmer_size };
    }
    /// otherwise, there is only one window
    else
    {
        _next_window = { _matrix, 0, 0 };
    }
}

impl::chained_window_iterator& impl::chained_window_iterator::operator++()
{
    _previous_window = std::move(_window);
    _window = std::move(_next_window);
    _next_window = _get_next_window();
    return *this;
}

window impl::chained_window_iterator::_get_next_window()
{
    /// Size of prefixes saved from the last window
    /// Only even k
    const auto prefix_size = _kmer_size / 2;
    const auto suffix_size = _kmer_size - prefix_size;

    /// continue the chain if possible
    if (size_t(_window.get_position() + suffix_size + _kmer_size) <= _matrix->width())
    {
        return { _matrix, _window.get_position() + suffix_size, _kmer_size };
    }
    /// if the chain is over, start the next one if possible
    else if ((_chain_start + 1 <= _last_chain_pos) && (_chain_start + 1 + _kmer_size <= _matrix->width()))
    {
        ++_chain_start;
        return { _matrix, _chain_start, _kmer_size };
    }
    /// otherwise, the iterator is over
    else
    {
        _kmer_size = 0;
        return { _matrix, 0, 0 };
    }
}


bool impl::chained_window_iterator::operator==(const chained_window_iterator& rhs) const noexcept
{
    return _window == rhs._window;
}

bool impl::chained_window_iterator::operator!=(const chained_window_iterator& rhs) const noexcept
{
    return !(*this == rhs);
}

std::tuple<window&, window&, window&> impl::chained_window_iterator::operator*() noexcept
{
    return { _previous_window, _window, _next_window };
}

to_windows::to_windows(const matrix* matrix, size_t kmer_size)
    : _matrix{ matrix }, _kmer_size{ kmer_size }//, _start_pos{ 0 }
{}

to_windows::const_iterator to_windows::begin() const
{
    return { _matrix, _kmer_size };
}

to_windows::const_iterator to_windows::end() const noexcept
{
    return { _matrix, 0 };
}


chain_windows::chain_windows(const matrix* matrix, size_t kmer_size)
    : _matrix{ matrix }, _kmer_size{ kmer_size }
{}

chain_windows::const_iterator chain_windows::begin() const
{
    return { _matrix, _kmer_size };
}

chain_windows::const_iterator chain_windows::end() const noexcept
{
    return { _matrix, 0 };
}