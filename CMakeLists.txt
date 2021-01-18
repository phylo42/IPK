cmake_minimum_required(VERSION 3.10)
project(xpas LANGUAGES CXX VERSION 0.2.0)

# default build type Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# define a hash map
if(NOT DEFINED HASH_MAP)
    set(HASH_MAP "USE_TSL_ROBIN_MAP")
endif()

if(${HASH_MAP} STREQUAL USE_SKA_FLAT_HASH_MAP)
    add_subdirectory(skarupke)
    set(HASH_MAP_EXPORT skarupke)
elseif(${HASH_MAP} STREQUAL USE_SKA_BYTELL_HASH_MAP)
    add_subdirectory(skarupke)
elseif(${HASH_MAP} STREQUAL USE_ABSL_FLAT_HASH_MAP)
    set(BUILD_TESTING OFF)
    add_subdirectory(abseil-cpp)
elseif(${HASH_MAP} STREQUAL USE_FOLLY_F14_FAST_MAP)
    add_subdirectory(folly)
elseif(${HASH_MAP} STREQUAL USE_PHMAP_FLAT_HASH_MAP)
    add_subdirectory(parallel-hashmap)
elseif(${HASH_MAP} STREQUAL USE_TSL_ROBIN_MAP)
    add_subdirectory(robin-map)
elseif(${HASH_MAP} STREQUAL USE_TSL_HOPSCOTCH_MAP)
    add_subdirectory(hopscotch-map)
endif()


# Catch2 includes
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/catch2/contrib")

add_subdirectory(catch2)
add_subdirectory(xpas)
add_subdirectory(build)
add_subdirectory(examples)
add_subdirectory(tests)
#add_subdirectory(pyxpas)