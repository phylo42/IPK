#include <algorithm>
#include "pk_compute.h"

using namespace xpas;
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

DCLA::DCLA(const window& window, size_t k, phylo_kmer::score_type eps)
    : _window(window)
    , _k(k)
    , _eps(eps)
{
    /// kmer_size can also be zero, which means the end() iterator
    const auto halfsize = size_t{ k / 2 };
    _prefix_size = (halfsize >= 1) ? halfsize : k;

    // preprocessing O(k): range product query
    preprocess();

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

        phylo_kmer::score_type eps_l = eps - best_score(j + h / 2, h - h / 2);
        phylo_kmer::score_type eps_r = eps - best_score(j, h / 2);

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

            //for (const auto& [a, a_score] : max)
            //{
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
                    if (score <= eps)
                    {
                        break;
                    }

                    phylo_kmer::key_type kmer;
                    if (prefix_sort)
                    {
                        kmer = (b << ((h - h / 2) * 2)) | a;
                    }
                    else
                    {
                        kmer = (a << ((h - h / 2) * 2)) | b;
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

void DCLA::preprocess()
{
    _best_scores = std::vector<phylo_kmer::score_type>(_k + 1, 1.0f);
    phylo_kmer::score_type product = 1.0f;
    for (size_t j = 0; j < _k; ++j)
    {
        const auto& [index_best, score_best] = _window.max_at(j);
        product += score_best;
        _best_scores[j + 1] = product;
    }

    /// We ignore 0s because they can not be highest values in the column
}

phylo_kmer::score_type DCLA::best_score(size_t start_pos, size_t h)
{
    return _best_scores[start_pos + h] - _best_scores[start_pos];
}


const std::vector<phylo_kmer>& DCLA::get_result() const
{
    return _result_list;
}