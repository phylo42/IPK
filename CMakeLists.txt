cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/catch2/contrib")

# exclude abseil-cpp tests
set(BUILD_TESTING OFF)

# hash maps
add_subdirectory(abseil-cpp)
add_subdirectory(folly)
add_subdirectory(skarupke)
add_subdirectory(parallel-hashmap)
add_subdirectory(robin-map)
add_subdirectory(hopscotch-map)
add_subdirectory(martinus)

add_subdirectory(catch2)
add_subdirectory(core)
add_subdirectory(examples)
add_subdirectory(tests)
add_subdirectory(utils)