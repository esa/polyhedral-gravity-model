# Exclude the YAML files from the python interface
set(PYTHON_INTERFACE_SRC ${SRC})

file(GLOB_RECURSE PYTHON_INTERFACE_SRC
    "${PROJECT_SOURCE_DIR}/src/polyhedralGravity/*.cpp"
    "${PROJECT_SOURCE_DIR}/src/polyhedralGravity/*.h"
    )
list(FILTER PYTHON_INTERFACE_SRC EXCLUDE REGEX ".*YAML.*")

# The OpenCL backend is only compiled in when it is enabled
if (NOT POLYHEDRAL_GRAVITY_ENABLE_OPENCL)
    list(FILTER PYTHON_INTERFACE_SRC EXCLUDE REGEX "/opencl/")
endif ()

pybind11_add_module(_core ${PYTHON_INTERFACE_SRC} PolyhedralGravityPython.cpp)

if (POLYHEDRAL_GRAVITY_ENABLE_OPENCL)
    target_link_libraries(_core PUBLIC polyhedral_gravity_opencl)
    # The generated kernel headers have to exist before the OpenCL sources are compiled
    add_dependencies(_core polyhedral_gravity_opencl_kernels)
endif ()

# Checks if MinGW is used, in this case, the .dll of the the standard library need to be linked statically to the
# Python interface since otherwise they won't be found
# https://github.com/pybind/pybind11/issues/2010 and https://cython.readthedocs.io/en/latest/src/tutorial/appendix.html
if(WIN32 AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    message(STATUS "Compiler is GNU, Platform is Windows, Linking libaries statically to the Python interface")
    target_link_libraries(_core PUBLIC
        -static-libgcc
        -static-libstdc++
        -W
        -Bstatic
        --whole-file
        -lwinpthread
        --no-whole-file
        )
endif()
