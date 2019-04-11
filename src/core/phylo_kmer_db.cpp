#include "core/phylo_kmer_db.h"

void phylo_kmer_db::put(kmer_t key, branch_node_t branch, score_t score)
{
    if (auto it = _map.find(key); it != _map.end())
    {
        if (auto inner_it = it->second.find(branch); inner_it != it->second.end())
        {
            if (inner_it->second < score)
            {
                _map[key][branch] = score;
            }
        }
        else
        {
            _map[key][branch] = score;
        }
    }
    else
    {
        _map[key][branch] = score;
    }
}


phylo_kmer_db::const_iterator phylo_kmer_db::begin() const
{
    return std::begin(_map);
}

phylo_kmer_db::const_iterator phylo_kmer_db::end() const
{
    return std::end(_map);
}

size_t phylo_kmer_db::size() const
{
    return _map.size();
}

