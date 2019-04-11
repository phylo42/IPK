cmake_minimum_required(VERSION 3.10)
project(rappas-core)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Werror")

set(SOURCE_FILES
        src/core/phylo_kmer.cpp
        src/core/phylo_kmer_db.cpp
)

include_directories(include)
include_directories(third-party/skarupke)

add_library(rappas-core ${SOURCE_FILES})
target_include_directories(rappas-core PUBLIC include)