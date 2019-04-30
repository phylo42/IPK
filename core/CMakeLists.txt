cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(rappas_core LANGUAGES CXX VERSION 0.1.3)

# WARNING: some dependencies are not listed here. See the top-level CMake file
# for more details

##########################################################################################
# PROTEINS
# WARNING: SEQ_TYPE_AA is not supported yet. We declare a target for this,
# but it can not be compiled. We also define this target EXCLUDE_FROM_ALL for a while
# to not build this target by default
add_library(core_p STATIC EXCLUDE_FROM_ALL "")
add_library(rappas::core_p ALIAS core_p)

target_compile_definitions(core_p
        PUBLIC
            SEQ_TYPE_AA)

target_sources(core_p
        PRIVATE
            src/core/phylo_kmer.cpp
            src/core/phylo_kmer_db.cpp
            src/core/kmer_iterator.cpp
            src/core/phylo_tree.cpp
        )

find_package(Boost REQUIRED COMPONENTS serialization)

target_link_libraries(core_p
        PUBLIC
            skarupke::flat_hash_map
            Boost::serialization
        PRIVATE
            rappas::utils
        )

target_include_directories(core_p
        PUBLIC
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Turn on the warnings and treat them as errors
target_compile_options(core_p
        PRIVATE
            -Wall -Wextra -Wpedantic -Werror
        )

# explicitly require c++17 standart
# TODO: check this out if it is possible to relax this to PUBLIC cxx_std_14 and PRIVATE cxx_std_17
target_compile_features(core_p
        PUBLIC
            cxx_std_17
        )


##########################################################################################
# DNA
# Create a library target and set properties
add_library(core_n STATIC "")
add_library(rappas::core_n ALIAS core_n)

target_compile_definitions(core_n
        PUBLIC
            SEQ_TYPE_DNA)

target_sources(core_n
        PRIVATE
            src/core/phylo_kmer.cpp
            src/core/phylo_kmer_db.cpp
            src/core/kmer_iterator.cpp
            src/core/phylo_tree.cpp
        )

find_package(Boost REQUIRED COMPONENTS serialization unit_test_framework)

target_link_libraries(core_n
        PUBLIC
            rappas::utils
            Boost::serialization
            Catch2::Catch2

            # hash maps
            skarupke::flat_hash_map
            absl::flat_hash_map
            folly_base
            phmap
            tsl::robin_map
            tsl::hopscotch_map
            martinus::robin_hood_map
        )

target_include_directories(core_n
        PUBLIC
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Turn on the warnings and treat them as errors
target_compile_options(core_n
        PRIVATE
            -Wall -Wextra -Werror #-Wpedantic
        )

# explicitly require c++17 standart
# TODO: check this out if it is possible to relax this to PUBLIC cxx_std_14 and PRIVATE cxx_std_17
target_compile_features(core_n
        PUBLIC
            cxx_std_17
        )

######################################################################################################
# Library installation
include(GNUInstallDirs)

# Exporting the core target
# 'make install' to the correct locations (provided by GNUInstallDirs)
# TODO: check if it works with 'make install'
install(TARGETS core_n
        EXPORT coreNTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(TARGETS core_n
       FILE ${CMAKE_BINARY_DIR}/rappasConfig.cmake)

export(PACKAGE rappas)