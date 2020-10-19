# xpas, the phylo k-mer database construction library


This is a core library used in RAPPAS2 and SHERPAS.

## Dependencies
xpas depends on:
- boost v1.67+ must be installed in your system
- third-party libraries included as submodules in the repository

## Install

1. Make sure boost v1.67+ is installed in your system.
2. Clone this repository with submodules.
3. Build it with one of hashmaps enabled (see below).

Example:
```
mkdir bin; cd bin
cmake -DHASH_MAP="USE_TSL_ROBIN_MAP -DCMAKE_CXX_FLAGS="-O3" ..
make
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
This library is not intended to be used stand-alone. It is used in RAPPAS2 and SHERPAS. Check out our [examples](https://github.com/phylo42/xpas/tree/master/examples) for the use cases.

