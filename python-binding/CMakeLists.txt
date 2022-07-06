# Include pybind11
include(pybind11)

# Stores the file of the python interface in a variable
set(pythonInterface ${CMAKE_CURRENT_SOURCE_DIR}/python-interface.cpp)

file(GLOB_RECURSE PYTHON_INTERFACE_SRC
        "${PROJECT_SOURCE_DIR}/src/polyhedralGravity/*.h"
        "${PROJECT_SOURCE_DIR}/src/polyhedralGravity/*.cpp")

list(FILTER PYTHON_INTERFACE_SRC EXCLUDE REGEX ".*YAML.*")

pybind11_add_module(polyhedral_gravity ${pythonInterface} ${PYTHON_INTERFACE_SRC})

target_include_directories(polyhedral_gravity PUBLIC
        "${PROJECT_SOURCE_DIR}/src"
        )

target_link_libraries(polyhedral_gravity PUBLIC
        spdlog::spdlog
        tetgen
        xsimd
        Thrust
        )

# Checks if MinGW is used, in this case, the .dll of the the standard library need to be linked statically to the
# Python interface since otherwise they won't be found
# https://github.com/pybind/pybind11/issues/2010 and https://cython.readthedocs.io/en/latest/src/tutorial/appendix.html
if(WIN32 AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    message(STATUS "Compiler is GNU, Platform is Windows, Linking libaries statically to the Python interface")
    target_link_libraries(polyhedral_gravity PUBLIC
        -static-libgcc
        -static-libstdc++
        -W
        -Bstatic
        --whole-file
        -lwinpthread
        --no-whole-file
        )
endif()