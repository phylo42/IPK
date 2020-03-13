cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

find_package(PythonInterp 3 REQUIRED)

message(STATUS "Python interpreter found: " ${PYTHON_VERSION_MAJOR})

find_package(Boost COMPONENTS python${PYTHON_VERSION_MAJOR})
find_package(PythonLibs 3 REQUIRED)
message(STATUS "Include dirs of Python: " ${PYTHON_INCLUDE_DIRS} )
message(STATUS "Libs of Python: " ${PYTHON_LIBRARIES} )


python_add_module(xpas)

target_sources(xpas PRIVATE src/pyxpas.cpp)

#file(COPY src/xpas.py DESTINATION .)
#add_test(NAME pyxpas COMMAND ${PYTHON_EXECUTABLE} xpas.py)


target_include_directories(xpas PUBLIC
        ${PYTHON_INCLUDE_DIRS}
)

target_link_libraries(xpas PRIVATE
        rappas::core_n
        Boost::python${PYTHON_VERSION_MAJOR}
        ${PYTHON_LIBRARIES}
)