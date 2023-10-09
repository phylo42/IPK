cmake_minimum_required(VERSION 3.10)
project(IPK LANGUAGES CXX VERSION 0.4.0)

# default build type Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Catch2 includes
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/catch2/contrib")

add_subdirectory(strasser)
add_subdirectory(indicators)
add_subdirectory(i2l)
add_subdirectory(ipk)
add_subdirectory(tools)

install(FILES ipk.py
        DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE)