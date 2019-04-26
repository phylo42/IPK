cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/catch2/contrib")

add_subdirectory(catch2)
add_subdirectory(core)
add_subdirectory(skarupke)
add_subdirectory(examples)
add_subdirectory(tests)