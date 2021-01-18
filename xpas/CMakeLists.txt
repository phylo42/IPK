cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

if(${HASH_MAP} STREQUAL "USE_SKA_FLAT_HASH_MAP")
    set(HASH_MAP_IMPLEMENTATION skarupke::flat_hash_map)
    set(HASH_MAP_EXPORT flat_hash_map)
elseif(${HASH_MAP} STREQUAL USE_SKA_BYTELL_HASH_MAP)
    message(WARNING "Using ska::flat_hash_map instead")
    set(HASH_MAP_IMPLEMENTATION skarupke::flat_hash_map)
    set(HASH_MAP_EXPORT flat_hash_map)
elseif(${HASH_MAP} STREQUAL USE_ABSL_FLAT_HASH_MAP)
    set(HASH_MAP_IMPLEMENTATION absl::flat_hash_map)
    set(HASH_MAP_EXPORT absl::flat_hash_map)
elseif(${HASH_MAP} STREQUAL USE_FOLLY_F14_FAST_MAP)
    set(HASH_MAP_IMPLEMENTATION folly_base)
    set(HASH_MAP_EXPORT folly_base)
elseif(${HASH_MAP} STREQUAL USE_PHMAP_FLAT_HASH_MAP)
    set(HASH_MAP_IMPLEMENTATION phmap)
    set(HASH_MAP_EXPORT phmap)
elseif(${HASH_MAP} STREQUAL USE_TSL_ROBIN_MAP)
    set(HASH_MAP_IMPLEMENTATION tsl::robin_map)
    set(HASH_MAP_EXPORT robin_map)
elseif(${HASH_MAP} STREQUAL USE_TSL_HOPSCOTCH_MAP)
    set(HASH_MAP_IMPLEMENTATION tsl::hopscotch_map)
    set(HASH_MAP_EXPORT hopscotch_map)
endif()


find_package(Boost REQUIRED COMPONENTS serialization iostreams)

set(SOURCES
        src/xpas/fasta.cpp include/xpas/fasta.h
        src/xpas/file_io.cpp include/xpas/file_io.h
        src/xpas/kmer_iterator.cpp include/xpas/kmer_iterator.h
        src/xpas/newick.cpp include/xpas/newick.h
        include/xpas/optional.h
        src/xpas/phylo_kmer.cpp include/xpas/phylo_kmer.h
        include/xpas/phylo_kmer_db.h
        src/xpas/phylo_node.cpp include/xpas/phylo_node.h
        src/xpas/phylo_tree.cpp include/xpas/phylo_tree.h
        src/xpas/seq.cpp include/xpas/seq.h
        src/xpas/seq_record.cpp include/xpas/seq_record.h
        include/xpas/serialization.h
        include/xpas/version.h)


set(LINK_LIBRARIES
        Boost::serialization Boost::iostreams
        ${HASH_MAP_IMPLEMENTATION})

set(INCLUDE_DIRECTORIES $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Turn on the warnings and treat them as errors
set(COMPILE_OPTIONS -Wall -Wextra -Wpedantic)

# explicitly require c++17 standard
# TODO: check this out if it is possible to relax this to PUBLIC cxx_std_14 and PRIVATE cxx_std_17
set(COMPILE_FEATURES cxx_std_17)

##########################################################################################
# Proteins
add_library(xpas_aa STATIC "")
add_library(xpas::aa ALIAS xpas_aa)

target_sources(xpas_aa PRIVATE ${SOURCES})
target_link_libraries(xpas_aa PUBLIC ${LINK_LIBRARIES})
target_include_directories(xpas_aa PUBLIC ${INCLUDE_DIRECTORIES})
target_compile_options(xpas_aa PRIVATE ${COMPILE_OPTIONS})
target_compile_features(xpas_aa PUBLIC ${COMPILE_FEATURES})
target_compile_definitions(xpas_aa
        PUBLIC
        ${HASH_MAP}
        SEQ_TYPE_AA)

##########################################################################################
# Proteins, positions included
add_library(xpas_aa_pos STATIC "")
add_library(xpas::aa_pos ALIAS xpas_aa_pos)

target_sources(xpas_aa_pos PRIVATE ${SOURCES})
target_link_libraries(xpas_aa_pos PUBLIC ${LINK_LIBRARIES})
target_include_directories(xpas_aa_pos PUBLIC ${INCLUDE_DIRECTORIES})
target_compile_options(xpas_aa_pos PRIVATE ${COMPILE_OPTIONS})
target_compile_features(xpas_aa_pos PUBLIC ${COMPILE_FEATURES})
target_compile_definitions(xpas_aa_pos
        PUBLIC
            ${HASH_MAP}
            SEQ_TYPE_AA
            KEEP_POSITIONS)

##########################################################################################
# DNA
# Create a library target and set properties
add_library(xpas_dna STATIC "")
add_library(xpas::dna ALIAS xpas_dna)

target_sources(xpas_dna PRIVATE ${SOURCES})
target_link_libraries(xpas_dna PUBLIC ${LINK_LIBRARIES})
target_include_directories(xpas_dna PUBLIC ${INCLUDE_DIRECTORIES})
target_compile_options(xpas_dna PRIVATE ${COMPILE_OPTIONS})
target_compile_features(xpas_dna PUBLIC ${COMPILE_FEATURES})
target_compile_definitions(xpas_dna
        PUBLIC
            ${HASH_MAP}
            SEQ_TYPE_DNA)

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

export(TARGETS xpas_dna ${HASH_MAP_EXPORT}
       FILE ${CMAKE_BINARY_DIR}/xpasConfig.cmake)

install(TARGETS xpas_aa
        EXPORT xpasAaTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(TARGETS xpas_aa ${HASH_MAP_EXPORT}
        FILE ${CMAKE_BINARY_DIR}/xpasConfig.cmake)

export(PACKAGE xpas)
