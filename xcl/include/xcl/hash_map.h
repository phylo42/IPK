#ifndef XPAS_HASH_MAP_H
#define XPAS_HASH_MAP_H

#ifdef USE_SKA_FLAT_HASH_MAP
#include <flat_hash_map/flat_hash_map.hpp>
#elif USE_SKA_BYTELL_HASH_MAP
#include <flat_hash_map/bytell_hash_map.hpp>
#elif USE_ABSL_FLAT_HASH_MAP
#include <absl/container/flat_hash_map.h>
#elif USE_FOLLY_F14_FAST_MAP
#include <folly/container/F14Map.h>
#include <folly/hash/Hash.h>
#elif USE_PHMAP_FLAT_HASH_MAP
#include <parallel_hashmap/phmap.h>
#elif USE_TSL_ROBIN_MAP
#include <tsl/robin_map.h>
#elif USE_TSL_HOPSCOTCH_MAP
#include <tsl/hopscotch_map.h>
#endif

namespace xpas
{
    template<typename... Args>
#ifdef USE_SKA_FLAT_HASH_MAP
    using hash_map = ska::flat_hash_map<Args...>;
#elif USE_SKA_BYTELL_HASH_MAP
    using hash_map = ska::bytell_hash_map<Args...>;
#elif USE_ABSL_FLAT_HASH_MAP
    using hash_map = absl::flat_hash_map<Args...>;
#elif USE_FOLLY_F14_FAST_MAP
    using hash_map = folly::F14FastMap<Args...>;
#elif USE_PHMAP_FLAT_HASH_MAP
    using hash_map = phmap::flat_hash_map<Args...>;
#elif USE_TSL_ROBIN_MAP
    using hash_map = tsl::robin_map<Args...>;
#elif USE_TSL_HOPSCOTCH_MAP
    using hash_map = tsl::hopscotch_map<Args...>;
#endif
/// requires manual intervention (does not support noexcept destructor and can not be defined this way)
///using hash_map = robin_hood::unordered_flat_map<Args...>;
///using hash_map = robin_hood::unordered_map<Args...>;
}
#endif //XPAS_HASH_MAP_H
