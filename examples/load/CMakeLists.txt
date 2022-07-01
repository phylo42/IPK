cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# Create an executable target and set properties
add_executable(xpas_example_load EXCLUDE_FROM_ALL "")

target_sources(xpas_example_load
    PRIVATE
        src/main.cpp
)

find_package(Boost REQUIRED COMPONENTS filesystem)

target_link_libraries(xpas_example_load
    PRIVATE
        xcl::dna
        Boost::filesystem
)

# Turn on the warnings and treat them as errors
target_compile_options(xpas_example_load
    PRIVATE
        -Wall -Wextra -Wpedantic -Werror
)

target_compile_features(xpas_example_load
    PRIVATE
        cxx_std_17
)