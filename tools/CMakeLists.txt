cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# Create an executable target and set properties
add_executable(diff-dna EXCLUDE_FROM_ALL "")
set_target_properties(diff-dna PROPERTIES OUTPUT_NAME xpas-diff-dna)
target_sources(diff-dna PRIVATE src/diff.cpp)
target_link_libraries(diff-dna PRIVATE xcl::dna)
target_compile_options(diff-dna PRIVATE -Wall -Wextra -Werror -Wpedantic)
set_property(TARGET diff-dna PROPERTY CXX_STANDARD 17)
target_compile_features(diff-dna PUBLIC cxx_std_17)


add_executable(diff-aa EXCLUDE_FROM_ALL "")
set_target_properties(diff-aa PROPERTIES OUTPUT_NAME xpas-diff-aa)
target_sources(diff-aa PRIVATE src/diff.cpp)
target_link_libraries(diff-aa PRIVATE xcl::aa)
target_compile_options(diff-aa PRIVATE -Wall -Wextra -Werror -Wpedantic)
set_property(TARGET diff-aa PROPERTY CXX_STANDARD 17)
target_compile_features(diff-aa PUBLIC cxx_std_17)
