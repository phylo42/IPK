cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# Create an executable target and set properties
add_executable(xpas_example_search EXCLUDE_FROM_ALL "")

target_sources(xpas_example_search
    PRIVATE
        src/main.cpp
)

target_link_libraries(xpas_example_search
    PRIVATE
        xpas::dna
)

# Turn on the warnings and treat them as errors
target_compile_options(xpas_example_search
    PRIVATE
        -Wall -Wextra -Wpedantic -Werror
)

target_compile_features(xpas_example_search
    PRIVATE
        cxx_std_17
)