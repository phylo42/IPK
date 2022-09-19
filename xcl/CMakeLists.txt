cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

if(${HASH_MAP} STREQUAL USE_SKA_FLAT_HASH_MAP)
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


find_package(Boost REQUIRED COMPONENTS serialization iostreams system filesystem)
find_package(ZLIB)

set(SOURCES
        src/xcl/fasta.cpp include/xcl/fasta.h
        src/xcl/file_io.cpp include/xcl/file_io.h
        include/xcl/hash_map.h
        src/xcl/kmer_iterator.cpp include/xcl/kmer_iterator.h
        src/xcl/newick.cpp include/xcl/newick.h
        include/xcl/optional.h
        src/xcl/phylo_kmer.cpp include/xcl/phylo_kmer.h
        include/xcl/phylo_kmer_db.h
        src/xcl/phylo_node.cpp include/xcl/phylo_node.h
        src/xcl/phylo_tree.cpp include/xcl/phylo_tree.h
        src/xcl/seq.cpp include/xcl/seq.h
        src/xcl/seq_record.cpp include/xcl/seq_record.h
        include/xcl/serialization.h
        include/xcl/version.h)


set(LINK_LIBRARIES
        Boost::serialization Boost::iostreams Boost::system Boost::filesystem
        ZLIB::ZLIB
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
add_library(xcl_aa STATIC "")
add_library(xcl::aa ALIAS xcl_aa)

target_sources(xcl_aa PRIVATE ${SOURCES})
target_link_libraries(xcl_aa PUBLIC ${LINK_LIBRARIES})
target_include_directories(xcl_aa PUBLIC ${INCLUDE_DIRECTORIES})
target_compile_options(xcl_aa PRIVATE ${COMPILE_OPTIONS})
target_compile_features(xcl_aa PUBLIC ${COMPILE_FEATURES})
target_compile_definitions(xcl_aa
        PUBLIC
        ${HASH_MAP}
        SEQ_TYPE_AA)

###########################################################################################
## Proteins, positions included
#add_library(xcl_aa_pos STATIC "")
#add_library(xcl::aa_pos ALIAS xcl_aa_pos)
#
#target_sources(xcl_aa_pos PRIVATE ${SOURCES})
#target_link_libraries(xcl_aa_pos PUBLIC ${LINK_LIBRARIES})
#target_include_directories(xcl_aa_pos PUBLIC ${INCLUDE_DIRECTORIES})
#set_property(TARGET xcl_aa_pos PROPERTY CXX_STANDARD 17)
#target_compile_options(xcl_aa_pos PRIVATE ${COMPILE_OPTIONS})
#target_compile_features(xcl_aa_pos PUBLIC ${COMPILE_FEATURES})
#target_compile_definitions(xcl_aa_pos
#        PUBLIC
#            ${HASH_MAP}
#            SEQ_TYPE_AA
#            KEEP_POSITIONS)

##########################################################################################
# DNA
# Create a library target and set properties
add_library(xcl_dna STATIC "")
add_library(xcl::dna ALIAS xcl_dna)

target_sources(xcl_dna PRIVATE ${SOURCES})
target_link_libraries(xcl_dna PUBLIC ${LINK_LIBRARIES})
target_include_directories(xcl_dna PUBLIC ${INCLUDE_DIRECTORIES})
target_compile_options(xcl_dna PRIVATE ${COMPILE_OPTIONS})
target_compile_features(xcl_dna PUBLIC ${COMPILE_FEATURES})
set_property(TARGET xcl_dna PROPERTY CXX_STANDARD 17)
target_compile_definitions(xcl_dna
        PUBLIC
            ${HASH_MAP}
            SEQ_TYPE_DNA
            KEEP_POSITIONS)


add_library(xcl_dna_pos STATIC "")
add_library(xcl::dna_pos ALIAS xcl_dna)

target_sources(xcl_dna_pos PRIVATE ${SOURCES})
target_link_libraries(xcl_dna_pos PUBLIC ${LINK_LIBRARIES})
target_include_directories(xcl_dna_pos PUBLIC ${INCLUDE_DIRECTORIES})
target_compile_options(xcl_dna_pos PRIVATE ${COMPILE_OPTIONS})
target_compile_features(xcl_dna_pos PUBLIC ${COMPILE_FEATURES})
set_property(TARGET xcl_dna_pos PROPERTY CXX_STANDARD 17)
target_compile_definitions(xcl_dna_pos
        PUBLIC
        ${HASH_MAP}
        SEQ_TYPE_DNA)

######################################################################################################
# Library installation
include(GNUInstallDirs)
#include(CMakePackageConfigHelpers)

# Exporting the core target
# 'make install' to the correct locations (provided by GNUInstallDirs)
# TODO: check if it works with 'make install'
install(TARGETS xcl_dna
        EXPORT xclDnaTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Give CMake access to the version of the library
#write_basic_package_version_file(
#        xclDnaConfigVersion.cmake
#        VERSION ${PACKAGE_VERSION}
#        COMPATIBILITY AnyNewerVersion
#)


#install(TARGETS xcl_aa
#        EXPORT xclAaTargets
#        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#export(TARGETS xcl_dna xcl_aa ${HASH_MAP_EXPORT}
export(TARGETS xcl_dna ${HASH_MAP_EXPORT}
        FILE ${CMAKE_BINARY_DIR}/xpasConfig.cmake)

export(PACKAGE xcl)
