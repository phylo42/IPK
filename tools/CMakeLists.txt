cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# Create an executable target and set properties
add_executable(diff EXCLUDE_FROM_ALL "")
set_target_properties(diff PROPERTIES OUTPUT_NAME xpas-diff)

target_sources(diff PRIVATE src/diff.cpp)

target_link_libraries(diff PRIVATE xpas::dna)

target_compile_options(diff PRIVATE -Wall -Wextra -Werror -Wpedantic)

set_property(TARGET diff PROPERTY CXX_STANDARD 17)

target_compile_features(diff PUBLIC cxx_std_17)
