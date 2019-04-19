cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# WARNING: some dependencies are not listed here. See the top-level CMake file
# for more details

# Create an executable target and set properties
add_executable(load_example EXCLUDE_FROM_ALL "")

target_sources(load_example
        PRIVATE
            src/main.cpp
        )

target_link_libraries(load_example
        PRIVATE
            rappas::core_n
        )

# Turn on the warnings and treat them as errors
target_compile_options(load_example
        PRIVATE
            -Wall -Wextra -Wpedantic -Werror
        )

target_compile_features(load_example
        PRIVATE
            cxx_std_17
        )