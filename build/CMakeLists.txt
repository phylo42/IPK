cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

######################################################################################################
# Dependencies
# WARNING: some dependencies are not listed here. See the top-level CMake file
# for more details
find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams)
find_package(OpenMP REQUIRED)

# The code and the libraries are the same for all targets.
# The only difference is the version of xpas we link against (xpas::dna, xpas::aa)
set(LINK_LIBRARIES Boost::program_options
        Boost::filesystem
        Boost::iostreams
        OpenMP::OpenMP_CXX
        strasser::csv_parser)


set(SOURCES
        src/alignment.cpp include/alignment.h
        src/ar.cpp include/ar.h
        src/branch_group.cpp include/branch_group.h
        src/command_line.cpp include/command_line.h
        src/db_builder.cpp include/db_builder.h
        src/exceptions.cpp include/exceptions.h
        src/extended_tree.cpp include/extended_tree.h
        src/filter.cpp include/filter.h
        src/main.cpp
        src/node_entry.cpp include/node_entry.h
        src/node_entry_view.cpp include/node_entry_view.h
        src/proba_matrix.cpp include/proba_matrix.h
        include/return.h
        include/row.h
        include/utils.h
        )

######################################################################################################
# Application target and properties
add_executable(build_dna "")
set_target_properties(build_dna PROPERTIES OUTPUT_NAME xpas-build-dna)

target_sources(build_dna PRIVATE ${SOURCES})

target_include_directories(build_dna
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
            ${CMAKE_CURRENT_SOURCE_DIR}/include/
            ${CMAKE_CURRENT_SOURCE_DIR}/strasser/
)

target_link_libraries(build_dna
        PRIVATE
            xpas::dna
            ${LINK_LIBRARIES}

        )

# Turn on the warnings and treat them as errors
target_compile_options(build_dna
    PRIVATE
        -Wall -Wextra #-Werror #-Wpedantic
    )

target_compile_features(build_dna
    PUBLIC
        cxx_std_17)

######################################################################################################
# Amino acids
add_executable(build_aa "")
set_target_properties(build_aa PROPERTIES OUTPUT_NAME xpas-build-aa)

target_sources(build_aa PRIVATE ${SOURCES})

target_include_directories(build_aa
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/
        ${CMAKE_CURRENT_SOURCE_DIR}/include/
        )

target_link_libraries(build_aa
        PRIVATE
        ${LINK_LIBRARIES}
        xpas::aa
        )

# Turn on the warnings and treat them as errors
target_compile_options(build_aa
        PRIVATE
        -Wall -Wextra #-Werror #-Wpedantic
        )

target_compile_features(build_aa
        PUBLIC
        cxx_std_17)

######################################################################################################
# Amino acids, positions included

add_executable(build_aa_pos "")
set_target_properties(build_aa_pos PROPERTIES OUTPUT_NAME xpas-build-aa-pos)

target_sources(build_aa_pos PRIVATE ${SOURCES})

target_include_directories(build_aa_pos
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/
        ${CMAKE_CURRENT_SOURCE_DIR}/include/
        )

target_link_libraries(build_aa_pos
        PRIVATE
        ${LINK_LIBRARIES}
        xpas::aa_pos
        )

# Turn on the warnings and treat them as errors
target_compile_options(build_aa_pos
        PRIVATE
        -Wall -Wextra #-Werror #-Wpedantic
        )

target_compile_features(build_aa_pos
        PUBLIC
        cxx_std_17)