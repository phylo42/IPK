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
        strasser::csv_parser
        indicators::indicators)


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
        src/window.cpp include/window.h
        src/pk_compute.cpp include/pk_compute.h
        src/proba_matrix.cpp include/proba_matrix.h
        include/return.h
        include/row.h
        include/utils.h
        )

######################################################################################################
# Application target and properties
add_executable(xpas_dna "")
set_target_properties(xpas_dna PROPERTIES OUTPUT_NAME xpas-dna)
target_sources(xpas_dna PRIVATE ${SOURCES})

target_include_directories(xpas_dna
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
            ${CMAKE_CURRENT_SOURCE_DIR}/include/
)

target_link_libraries(xpas_dna
        PUBLIC
            ${LINK_LIBRARIES}
            xcl::dna

        )

# Turn on the warnings and treat them as errors
target_compile_options(xpas_dna
    PRIVATE
        -Wall -Wextra #-Werror #-Wpedantic
    )


set_property(TARGET xpas_dna PROPERTY CXX_STANDARD 17)
target_compile_features(xpas_dna
    PUBLIC
        cxx_std_17)

######################################################################################################
## Amino acids
#add_executable(xpas_aa "")
#set_target_properties(xpas_aa PROPERTIES OUTPUT_NAME xpas-aa)
#target_sources(xpas_aa PRIVATE ${SOURCES})
#
#target_include_directories(xpas_aa
#        PRIVATE
#        ${CMAKE_CURRENT_SOURCE_DIR}/src/
#        ${CMAKE_CURRENT_SOURCE_DIR}/include/
#        )
#
#target_link_libraries(xpas_aa
#        PRIVATE
#        ${LINK_LIBRARIES}
#        xcl::aa
#        )
#
## Turn on the warnings and treat them as errors
#target_compile_options(xpas_aa
#        PRIVATE
#        -Wall -Wextra #-Werror #-Wpedantic
#        )
#
#target_compile_features(xpas_aa
#        PUBLIC
#        cxx_std_17)

######################################################################################################
# Amino acids, positions included

#add_executable(xpas_aa_pos "")
#
#target_sources(xpas_aa_pos PRIVATE ${SOURCES})
#
#target_include_directories(xpas_aa_pos
#        PRIVATE
#        ${CMAKE_CURRENT_SOURCE_DIR}/src/
#        ${CMAKE_CURRENT_SOURCE_DIR}/include/
#        )
#
#target_link_libraries(xpas_aa_pos
#        PRIVATE
#        ${LINK_LIBRARIES}
#        xcl::aa_pos
#        )
#
## Turn on the warnings and treat them as errors
#target_compile_options(xpas_aa_pos
#        PRIVATE
#        -Wall -Wextra #-Werror #-Wpedantic
#        )
#
#set_property(TARGET xpas_aa_pos PROPERTY CXX_STANDARD 17)
#target_compile_features(xpas_aa_pos
#        PUBLIC
#        cxx_std_17)