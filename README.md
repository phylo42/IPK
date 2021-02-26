# xpas, the phylo k-mer database construction library

xpas is the shared part of RAPPAS2 and SHERPAS, which builds phylo k-mer databases and provides efficient access to previously built databases.

## Dependencies
xpas depends on:
- boost v1.67+, boost-serialization, boost-filesystem, boost-iostreams, boost-program-options
- zlib
- third-party libraries included as submodules in the repository

## Install

1. Make sure boost libraries are installed in your system.
2. Clone this repository with submodules.
3. Build it with one of hashmaps enabled (see below).

Example:
```
# Install dependencies
sudo apt-get update && sudo apt-get install -yq libboost-dev libboost-serialization-dev libboost-filesystem-dev libboost-iostreams-dev libboost-program-options-dev zlib1g-dev

# Clone the repo
git clone --recursive https://github.com/phylo42/xpas.git

# Build xpas
cd xpas; mkdir bin; cd bin
cmake -DHASH_MAP=USE_TSL_ROBIN_MAP -DCMAKE_CXX_FLAGS="-O3" ..
make -j4
```

HASH_MAP options:
- USE_SKA_FLAT_HASH_MAP
- USE_SKA_BYTELL_HASH_MAP
- USE_ABSL_FLAT_HASH_MAP
- USE_FOLLY_F14_FAST_MAP
- USE_PHMAP_FLAT_HASH_MAP
- USE_TSL_ROBIN_MAP
- USE_TSL_HOPSCOTCH_MAP

## Usage

To build a database:
```
python xpas.py build -s [nucl|amino] -b `which phyml` -w workdir -r alignment.fasta -t tree.newick -k 10
```

To check out all options:
```
python xpas.py build --help
```

Databases can be used in [RAPPAS2](https://github.com/phylo42/rappas2) and [SHERPAS](https://github.com/phylo42/sherpas). Also check out our [examples](https://github.com/phylo42/xpas/tree/master/examples) to see how to use xpas to retrieve information from databases.
