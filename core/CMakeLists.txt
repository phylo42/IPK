cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# WARNING: some dependencies are not listed here. See the top-level CMake file
# for more details

# Create a library target and set properties
add_library(core STATIC "")
add_library(rappas::core ALIAS core)

target_sources(core
        PRIVATE
            src/core/phylo_kmer.cpp
            src/core/phylo_kmer_db.cpp
        )

# WARNING: third-party dependency here
target_link_libraries(core
        PUBLIC
            skarupke::flat_hash_map
        )

target_include_directories(core
        PUBLIC
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Turn on the warnings and treat them as errors
target_compile_options(core
        PRIVATE
            -Wall -Wextra -Wpedantic -Werror
        )

# explicitly require c++17 standart
# TODO: check this out if it is possible to relax this to PUBLIC cxx_std_14 and PRIVATE cxx_std_17
target_compile_features(core
        PUBLIC
            cxx_std_17
        )


######################################################################################################
# Library installation
include(GNUInstallDirs)

# Exporting the core target
# 'make install' to the correct locations (provided by GNUInstallDirs)
# TODO: check if it works with 'make install'
install(TARGETS core
        EXPORT coreTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(TARGETS core
        FILE ${CMAKE_BINARY_DIR}/rappasConfig.cmake)

export(PACKAGE rappas)