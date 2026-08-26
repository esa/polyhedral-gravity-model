# Including the actual Polyhedral source files
file(GLOB_RECURSE SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_library(${PROJECT_NAME}_lib OBJECT ${SRC} Info.h)

# Kokkos compiles the device code as CXX and injects the paradigm's flags (including
# -expt-relaxed-constexpr, see Kokkos_ENABLE_CUDA_CONSTEXPR in cmake/kokkos.cmake) through its
# interface target, so no device specific compile options have to be set here.
