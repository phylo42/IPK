cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(xpas LANGUAGES CXX VERSION 0.1.7)

if(${HASH_MAP} STREQUAL "USE_SKA_FLAT_HASH_MAP")
    set(HASH_MAP_IMPLEMENTATION skarupke::flat_hash_map)
elseif(${HASH_MAP} STREQUAL USE_SKA_BYTELL_HASH_MAP)
    set(HASH_MAP_IMPLEMENTATION skarupke::flat_hash_map)
elseif(${HASH_MAP} STREQUAL USE_ABSL_FLAT_HASH_MAP)
    set(HASH_MAP_IMPLEMENTATION absl::flat_hash_map)
elseif(${HASH_MAP} STREQUAL USE_FOLLY_F14_FAST_MAP)
    set(HASH_MAP_IMPLEMENTATION folly_base)
elseif(${HASH_MAP} STREQUAL USE_PHMAP_FLAT_HASH_MAP)
    set(HASH_MAP_IMPLEMENTATION phmap)
elseif(${HASH_MAP} STREQUAL USE_TSL_ROBIN_MAP)
    set(HASH_MAP_IMPLEMENTATION tsl::robin_map)
elseif(${HASH_MAP} STREQUAL USE_TSL_HOPSCOTCH_MAP)
    set(HASH_MAP_IMPLEMENTATION tsl::hopscotch_map)
endif()

##########################################################################################
# PROTEINS
# WARNING: SEQ_TYPE_AA is not supported yet. We declare a target for this,
# but it can not be compiled. We also define this target EXCLUDE_FROM_ALL for a while
# to not build this target by default
add_library(xpas_aa STATIC "")
add_library(xpas::aa ALIAS xpas_aa)

target_compile_definitions(xpas_aa
        PUBLIC
            ${HASH_MAP}
            SEQ_TYPE_AA)

target_sources(xpas_aa
        PRIVATE
        src/xpas/phylo_kmer.cpp include/xpas/phylo_kmer.h
        src/xpas/phylo_kmer_db.cpp include/xpas/phylo_kmer_db.h
        src/xpas/kmer_iterator.cpp include/xpas/kmer_iterator.h
        src/xpas/phylo_tree.cpp include/xpas/phylo_tree.h
        src/xpas/newick.cpp include/xpas/newick.h
        include/xpas/serialization.h
        include/xpas/version.h
        )

find_package(Boost REQUIRED COMPONENTS serialization unit_test_framework)

target_link_libraries(xpas_aa
        PUBLIC
        xpas::utils
        Boost::serialization
        ${HASH_MAP_IMPLEMENTATION}
        )

target_include_directories(xpas_aa
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Turn on the warnings and treat them as errors
target_compile_options(xpas_aa
        PRIVATE
        -Wall -Wextra -Werror #-Wpedantic
        )

# explicitly require c++17 standard
target_compile_features(xpas_aa
        PUBLIC
        cxx_std_17
        )

##########################################################################################
# DNA
# Create a library target and set properties
add_library(xpas_dna STATIC "")
add_library(xpas::dna ALIAS xpas_dna)

target_compile_definitions(xpas_dna
        PUBLIC
            ${HASH_MAP}
            SEQ_TYPE_DNA)

target_sources(xpas_dna
        PRIVATE
            src/xpas/phylo_kmer.cpp include/xpas/phylo_kmer.h
            src/xpas/phylo_kmer_db.cpp include/xpas/phylo_kmer_db.h
            src/xpas/kmer_iterator.cpp include/xpas/kmer_iterator.h
            src/xpas/phylo_tree.cpp include/xpas/phylo_tree.h
            src/xpas/newick.cpp include/xpas/newick.h
            include/xpas/serialization.h
            include/xpas/version.h
        )

find_package(Boost REQUIRED COMPONENTS serialization unit_test_framework)


target_link_libraries(xpas_dna
        PUBLIC
            xpas::utils
            Boost::serialization
            ${HASH_MAP_IMPLEMENTATION}
        )

target_include_directories(xpas_dna
        PUBLIC
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Turn on the warnings and treat them as errors
target_compile_options(xpas_dna
        PRIVATE
            -Wall -Wextra -Werror #-Wpedantic
        )

# explicitly require c++17 standard
# TODO: check this out if it is possible to relax this to PUBLIC cxx_std_14 and PRIVATE cxx_std_17
target_compile_features(xpas_dna
        PUBLIC
            cxx_std_17
        )

######################################################################################################
# Library installation
include(GNUInstallDirs)

# Exporting the core target
# 'make install' to the correct locations (provided by GNUInstallDirs)
# TODO: check if it works with 'make install'
install(TARGETS xpas_dna
        EXPORT xpasDnaTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(TARGETS xpas_dna
       FILE ${CMAKE_BINARY_DIR}/xpasConfig.cmake)

install(TARGETS xpas_aa
        EXPORT xpasAaTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(TARGETS xpas_aa
        FILE ${CMAKE_BINARY_DIR}/xpasConfig.cmake)

export(PACKAGE xpas)