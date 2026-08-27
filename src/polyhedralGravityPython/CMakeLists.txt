# Exclude the YAML files from the python interface
set(PYTHON_INTERFACE_SRC ${SRC})

file(GLOB_RECURSE PYTHON_INTERFACE_SRC
    "${PROJECT_SOURCE_DIR}/src/polyhedralGravity/*.cpp"
    "${PROJECT_SOURCE_DIR}/src/polyhedralGravity/*.h"
    )
list(FILTER PYTHON_INTERFACE_SRC EXCLUDE REGEX ".*YAML.*")

# pybind11 adds link time optimization to the module by default. That cannot be combined with a Kokkos
# device backend: every object then also carries the device code as a fatbin, and the LTO stage re-emits
# those into a single assembly file, where they collide with
# "Error: symbol `fatbinData' is already defined". NO_EXTRAS switches pybind11's LTO off; the stripping,
# which is the other half of the extras, is requested explicitly below.
if (Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_HIP OR Kokkos_ENABLE_SYCL)
    message(STATUS "Building the Python interface without link time optimization, as it is incompatible "
            "with the Kokkos device backend")
    pybind11_add_module(_core NO_EXTRAS ${PYTHON_INTERFACE_SRC} PolyhedralGravityPython.cpp)
    if (NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|RelWithDebInfo|None)$")
        pybind11_strip(_core)
    endif ()
else ()
    pybind11_add_module(_core ${PYTHON_INTERFACE_SRC} PolyhedralGravityPython.cpp)
endif ()

# scikit-build-core builds the wheel from what CMake installs, so the module has to be installed next to
# the pure Python part of the package (see wheel.packages in pyproject.toml). SKBUILD is only defined when
# the build was started by pip; a plain CMake build keeps the module in the build tree as before.
if (SKBUILD)
    # The COMPONENT is what keeps the wheel to just this module: the bundled dependencies (Kokkos above all)
    # carry install rules of their own, and a plain "cmake --install" would put their headers, static
    # libraries, and tools into site-packages. pyproject.toml installs this component only.
    install(TARGETS _core LIBRARY DESTINATION polyhedral_gravity COMPONENT PolyhedralGravityPython)
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
