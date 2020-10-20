cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# Tests
enable_testing()

add_executable(test_dna src/main.cpp)

find_package(Boost REQUIRED COMPONENTS filesystem iostreams)

target_link_libraries(test_dna
        PRIVATE
            xpas::dna
            Catch2::Catch2
            Boost::filesystem
            Boost::iostreams
        )

include(CTest)
include(Catch)
catch_discover_tests(test_dna)