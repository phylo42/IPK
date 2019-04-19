cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# WARNING: some dependencies are not listed here. See the top-level CMake file
# for more details

# Create an executable target and set properties
add_executable(basic_example EXCLUDE_FROM_ALL)

target_sources(basic_example
        PRIVATE
            src/main.cpp
        )

target_link_libraries(basic_example
        PRIVATE
            rappas::core_dna
        )

# Turn on the warnings and treat them as errors
target_compile_options(basic_example
        PRIVATE
            -Wall -Wextra -Wpedantic -Werror
        )

target_compile_features(basic_example
        PRIVATE
            cxx_std_17
        )