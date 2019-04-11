cmake_minimum_required(VERSION 3.13)
project(rappas_core)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Werror")

set(SOURCE_FILES
        src/phylo_kmer.cpp
        src/phylo_kmer_db.cpp
)

include_directories(include)
include_directories(third-party/skarupke)

add_library(rappas_core ${SOURCE_FILES})
target_include_directories(rappas_core PUBLIC include)