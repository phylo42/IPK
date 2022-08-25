#include <algorithm>
#include "pk_compute.h"

using namespace xpas;
using namespace xcl;
using xpas::impl::vector_type;

bool kmer_score_comparator(const phylo_kmer& k1, const phylo_kmer& k2)
{
    return k1.score > k2.score;
}

/// Creates a vector of 1-mers from a column of PP matrix
std::vector<phylo_kmer> as_column(const window& window, size_t j, phylo_kmer::score_type eps)
{
    std::vector<phylo_kmer> column;
    for (size_t i = 0; i < seq_traits::alphabet_size; ++i)
    {
        const auto& element = window.get(i, j);
        if (element > eps)
        {
            column.push_back({ static_cast<phylo_kmer::key_type>(i), element });
        }
    }
    return column;
}

DCLA::DCLA(const window& window, size_t k)
    : _window(window)
    , _k(k)
{
}


void DCLA::run(phylo_kmer::score_type eps)
{
    _result_list = DC(0, _k, eps);
}


// j is the starat position of the window
// h is the length of the window
std::vector<phylo_kmer> DCLA::DC(size_t j, size_t h, phylo_kmer::score_type eps)
{
    // trivial case
    if (h == 1)
    {
        return as_column(_window, j, eps);
    }
    else
    {
        std::vector<phylo_kmer> result_vector;
        std::vector<phylo_kmer>& result = (h == _k) ? _result_list : result_vector;

        phylo_kmer::score_type eps_l = eps - _window.range_max_product(j + h / 2, h - h / 2);
        phylo_kmer::score_type eps_r = eps - _window.range_max_product(j, h / 2);

        //phylo_kmer::score_type eps_l = eps / best_score(j + h / 2, h - h / 2);
        //phylo_kmer::score_type eps_r = eps / best_score(j, h / 2);

        auto l = DC(j, h / 2, eps_l);
        auto r = DC(j + h / 2, h - h / 2, eps_r);

        // let's sort not suffixes, but whichever is less to sort, suffixes or prefixes
        bool prefix_sort = l.size() < r.size();
        auto min = prefix_sort ? l : r;
        auto max = prefix_sort ? r : l;

        auto eps_min = prefix_sort ? eps_l : eps_r;
        auto eps_max = prefix_sort ? eps_r : eps_l;

        if (!min.empty())
        {
            std::sort(min.begin(), min.end(), kmer_score_comparator);

            size_t i = 0;
            while (i < max.size())
            {
                const auto& [a, a_score] = max[i];
                if (a_score < eps_max)
                {
                    break;
                }

                size_t i2 = 0;
                while (i2 < min.size())
                {
                    const auto& [b, b_score] = min[i2];
                    if (b_score < eps_min)
                    {
                        break;
                    }

                    const auto score = a_score + b_score;
                    //const auto score = a_score * b_score;
                    if (score <= eps)
                    {
                        break;
                    }

                    phylo_kmer::key_type kmer;
                    if (prefix_sort)
                    {
                        kmer = (b << ((h - h / 2) * bit_length<seq_type>())) | a;
                    }
                    else
                    {
                        kmer = (a << ((h - h / 2) * bit_length<seq_type>())) | b;
                    }
                    result.push_back({ kmer, score });

                    i2++;
                }
                i++;
            }
        }
        return result;
    }
}
/*
void DCLA::preprocess()
{
    _best_scores = std::vector<phylo_kmer::score_type>(_k + 1, 1.0f);
    phylo_kmer::score_type product = 1.0f;
    for (size_t j = 0; j < _k; ++j)
    {
        const auto& [index_best, score_best] = _window.max_at(j);
        product += score_best;
        //product *= score_best;
        _best_scores[j + 1] = product;
    }

    /// We ignore 0s because they can not be highest values in the column
}

phylo_kmer::score_type DCLA::best_score(size_t start_pos, size_t h)
{
    return _best_scores[start_pos + h] - _best_scores[start_pos];
    //return _best_scores[start_pos + h] / _best_scores[start_pos];
}*/


const std::vector<phylo_kmer>& DCLA::get_result() const
{
    return _result_list;
}

DCCW::DCCW(const window& window, std::vector<phylo_kmer>& prefixes,
           size_t k, phylo_kmer::score_type lookbehind, phylo_kmer::score_type lookahead)
    : _window(window)
    , _k(k)
    , _lookahead(lookahead)
    , _lookbehind(lookbehind)
    , _prefixes(prefixes)
{
}


void DCCW::run(phylo_kmer::score_type eps)
{
    DCLA dc(_window, _k);

    phylo_kmer::score_type eps_r = eps - _window.range_max_product(0, _k / 2);
    phylo_kmer::score_type eps_l = eps - _window.range_max_product(_k / 2, _k - _k / 2);

    auto& L = _prefixes;
    if (L.empty())
    {
        L = dc.DC(0, _k / 2, eps_l);
    }

    // Copy elision
    _suffixes = dc.DC(_k / 2, _k - _k / 2, std::min(eps_r, eps - _lookahead));
    auto& R = _suffixes;

    // Let's sort not suffixes, but whichever is less to sort, suffixes or prefixes
    // The trick is that L can contain more dead prefixes (which were alive suffixes in the previous window),
    // and R can contain dead suffixes (which will be alive prefixes in the next window).
    // Let's first find how many alive prefixes and suffixes we have in the current window.
    auto last_prefix = L.end();
    auto last_suffix = R.end();

    // The best prefix score of the previous window was better than the best suffix score of
    // the current window. Then, more strings of L are alive in W_prev than in W.
    // => Need to partition L to find only the part of strings that are alive in W.
    if (eps - _lookbehind < eps_l)
    {
        auto is_alive_prefix = [eps_l](const auto& pk) { return pk.score > eps_l; };

        last_prefix = std::partition(L.begin(), L.end(), is_alive_prefix);
    }

    // The same for strings of R and the lookahead score for the next window
    if (eps - _lookahead < eps_r)
    {
        auto is_alive_suffix = [eps_r](const auto& pk) { return pk.score > eps_r; };
        last_suffix = std::partition(R.begin(), R.end(), is_alive_suffix);
    }

    size_t num_alive_prefixes = std::distance(L.begin(), last_prefix);
    size_t num_alive_suffixes = std::distance(R.begin(), last_suffix);

    bool prefix_sort = num_alive_prefixes < num_alive_suffixes;
    auto& min = prefix_sort ? L : R;
    auto& max = prefix_sort ? R : L;
    auto last_min = prefix_sort ? last_prefix : last_suffix;

    if (!min.empty())
    {
        auto eps_min = prefix_sort ? eps_l : eps_r;
        auto eps_max = prefix_sort ? eps_r : eps_l;

        std::sort(min.begin(), last_min, kmer_score_comparator);

        size_t i = 0;
        while (i < max.size())
        {
            const auto& [a, a_score] = max[i];
            if (a_score < eps_max)
            {
                break;
            }

            size_t j = 0;
            while (j < min.size())
            {
                const auto& [b, b_score] = min[j];
                if (b_score < eps_min)
                {
                    break;
                }

                const auto score = a_score * b_score;
                if (score <= eps)
                {
                    break;
                }

                phylo_kmer::key_type kmer;
                if (prefix_sort)
                {
                    kmer = (b << ((_k - _k / 2) * 2)) | a;
                }
                else
                {
                    kmer = (a << ((_k - _k / 2) * 2)) | b;
                }
                _result_list.push_back({ kmer, score });

                j++;
            }

            i++;
        }
    }
}


const std::vector<phylo_kmer>& DCCW::get_result() const
{
    return _result_list;
}

std::vector<phylo_kmer>&& DCCW::get_suffixes()
{
    return std::move(_suffixes);
}



BB::BB(const window& window, size_t k, phylo_kmer::score_type eps)
    : _window(window)
      , _k(k)
      , _best_suffix_score()
{
    if (_window.empty())
    {
        throw std::runtime_error("The matrix is empty.");
    }

    if (_window.size() != k)
    {
        throw std::runtime_error("The size of the window is not k");
    }

    preprocess();

    run(eps);
}

void BB::run(phylo_kmer::score_type eps)
{
    // Recursive BB
    for (size_t i = 0; i < seq_traits::alphabet_size; ++i)
    {
        bb(i, 0, 0, 1.0, eps);
    }
}

void BB::bb(size_t i, size_t j, phylo_kmer::key_type prefix, phylo_kmer::score_type score, phylo_kmer::score_type eps)
{
    score = score + _window.get(i, j);
    //score = score * _window.get(i, j);
    prefix = (prefix << 2) | i;

    if (j == _k - 1)
    {
        if (score > eps)
        {
            _result_list.push_back({prefix, score});
            //_map[prefix] = score;
            return;
        }
        else
        {
            return;
        }
    }

    const auto best_suffix = _best_suffix_score[_k - (j + 2)];
    if (score + best_suffix <= eps)
    //if (score * best_suffix <= eps)
    {
        return;
    }
    else
    {
        for (size_t i2 = 0; i2 < seq_traits::alphabet_size; ++i2)
        {
            bb(i2, j + 1, prefix, score, eps);
        }
        return;
    }
}

const std::vector<phylo_kmer>& BB::get_result() const
{
    return _result_list;
}

void BB::preprocess()
{
    // precalc the scores of the best suffixes
    phylo_kmer::key_type prefix = 0;
    //phylo_kmer::score_type score = 1.0;
    phylo_kmer::score_type score = 0.0;

    for (size_t i = 0; i < _k; ++i)
    {
        const auto& [index_best, score_best] = _window.max_at(_k - i - 1);

        prefix = (index_best << i * 2) | prefix;
        //score = score * score_best;
        score = score + score_best;

        _best_suffix_score.push_back(score);
    }
}


