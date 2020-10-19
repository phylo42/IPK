cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# Tests
include(CTest)
include(Catch)
enable_testing()

add_executable(test_dna EXCLUDE_FROM_ALL "")

target_sources(test_dna
        PRIVATE
            src/main.cpp
        )

find_package(Boost REQUIRED COMPONENTS filesystem)

target_link_libraries(test_dna
        PRIVATE
            xpas::dna
            Catch2::Catch2
            Boost::filesystem
        )

catch_discover_tests(test_dna)