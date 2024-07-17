# IPK: Inference of Phylo-K-mers
[![build](https://github.com/phylo42/IPK/actions/workflows/build.yml/badge.svg)](https://github.com/phylo42/IPK/actions/workflows/build.yml)

**Please cite:**  [![doi](https://img.shields.io/static/v1?label=doi&message=10.1093/bioinformatics/btad692&color=blue)](https://doi.org/10.1093/bioinformatics/btad692) [1]

IPK is a tool for computing phylo-k-mers for a fixed phylogeny.

[Link to documentation](https://phylo-k-mers.readthedocs.io/) 

For more information about phylo-k-mers, see our papers: [one](https://academic.oup.com/bioinformatics/article/39/12/btad692/7425449),  [two](https://doi.org/10.1093/bioinformatics/btz068), [three](https://doi.org/10.1093/bioinformatics/btaa1020).

If you want to experiment with phylo-k-mers and write some code, check out our [examples](https://github.com/phylo42/I2L/tree/master/examples) to see how to use the code of IPK.


## Rapid installation via Bioconda
It is advised to install the package in a new environment, because our C++ dependencies are strict and may clash with other packages requiring, for instance, libboost.
We also recommend to use `mamba, which is more faster in solving environment dependencies.
```
conda create -n my_env
conda activate my_env
conda config --set channel_priority strict
mamba install ipk
```

## References

---
[1] Romashchenko, Nikolai et al. EPIK: precise and scalable evolutionary placement with informative k-mers. Bioinformatics, 39.12 (2023), btad692.

[2] Romashchenko, Nikolai. Computing informative k-mers for phylogenetic placement. Diss. Université Montpellier, 2021.

[3] Linard, Benjamin et al. "Rapid alignment-free phylogenetic identification of metagenomic sequences." Bioinformatics 35.18 (2019): 3303-3312.

[4] Scholz, Guillaume E., et al. "Rapid screening and detection of inter-type viral recombinants using phylo-k-mers." Bioinformatics 36.22-23 (2020): 5351-5360.
