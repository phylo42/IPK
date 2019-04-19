cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# WARNING: some dependencies are not listed here. See the top-level CMake file
# for more details

##########################################################################################
# PROTEINS
# WARNING: SEQ_TYPE_PROTEIN is not supported yet. We declare a target for this,
# but it will not compile. We also define this target EXCLUDE_FROM_ALL for a while
# to not build this target by default
add_library(core_aa STATIC EXCLUDE_FROM_ALL)
add_library(rappas::core_aa ALIAS core_aa)

target_compile_definitions(core_aa
        PUBLIC
            SEQ_TYPE_AA)

target_sources(core_aa
        PRIVATE
            src/core/phylo_kmer.cpp
            src/core/phylo_kmer_db.cpp
        )

find_package(Boost REQUIRED COMPONENTS serialization)

target_link_libraries(core_aa
        PUBLIC
            skarupke::flat_hash_map
            Boost::serialization
        )

target_include_directories(core_aa
        PUBLIC
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Turn on the warnings and treat them as errors
target_compile_options(core_aa
        PRIVATE
            -Wall -Wextra -Wpedantic -Werror
        )

# explicitly require c++17 standart
# TODO: check this out if it is possible to relax this to PUBLIC cxx_std_14 and PRIVATE cxx_std_17
target_compile_features(core_aa
        PUBLIC
            cxx_std_17
        )


##########################################################################################
# DNA
# Create a library target and set properties
add_library(core_dna STATIC "")
add_library(rappas::core_dna ALIAS core_dna)

target_compile_definitions(core_dna
        PUBLIC
            SEQ_TYPE_DNA)

target_sources(core_dna
        PRIVATE
            src/core/phylo_kmer.cpp
            src/core/phylo_kmer_db.cpp
        )

find_package(Boost REQUIRED COMPONENTS serialization)

target_link_libraries(core_dna
        PUBLIC
            skarupke::flat_hash_map
            Boost::serialization
        )

target_include_directories(core_dna
        PUBLIC
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Turn on the warnings and treat them as errors
target_compile_options(core_dna
        PRIVATE
            -Wall -Wextra -Wpedantic -Werror
        )

# explicitly require c++17 standart
# TODO: check this out if it is possible to relax this to PUBLIC cxx_std_14 and PRIVATE cxx_std_17
target_compile_features(core_dna
        PUBLIC
            cxx_std_17
        )


######################################################################################################
# Library installation
include(GNUInstallDirs)

# Exporting the core target
# 'make install' to the correct locations (provided by GNUInstallDirs)
# TODO: check if it works with 'make install'
install(TARGETS core_dna
        EXPORT coreDnaTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(TARGETS core_dna
        FILE ${CMAKE_BINARY_DIR}/rappasConfig.cmake)

export(PACKAGE rappas)