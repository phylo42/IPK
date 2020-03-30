cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

find_package(PythonInterp 3 REQUIRED)
find_package(Boost COMPONENTS python3 numpy3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")



python_add_module(pyxpas)

target_sources(pyxpas PRIVATE src/pyxpas.cpp)

#file(COPY src/xpas.py DESTINATION .)
#add_test(NAME pyxpas COMMAND ${PYTHON_EXECUTABLE} xpas.py)


target_include_directories(pyxpas PUBLIC
        ${PYTHON_INCLUDE_DIRS}
)

target_link_libraries(pyxpas PRIVATE
        rappas::core_n
        #Boost::boost
        Boost::python3
        Boost::numpy3
        ${PYTHON_LIBRARIES}
)