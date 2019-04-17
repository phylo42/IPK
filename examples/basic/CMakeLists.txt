cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# WARNING: some dependencies are not listed here. See the top-level CMake file
# for more details

# Create an executable target and set properties
add_executable(basic EXCLUDE_FROM_ALL)

target_sources(basic
        PRIVATE
            src/main.cpp
        )
target_include_directories(basic
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
        )

target_link_libraries(basic
        PRIVATE
            rappas::core
        )

# Turn on the warnings and treat them as errors
target_compile_options(basic
        PRIVATE
            -Wall -Wextra -Wpedantic -Werror
        )

target_compile_features(basic
        PRIVATE
            cxx_std_17
        )