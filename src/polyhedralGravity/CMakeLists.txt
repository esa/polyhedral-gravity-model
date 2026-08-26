# Including the actual Polyhedral source files
file(GLOB_RECURSE SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_library(${PROJECT_NAME}_lib OBJECT ${SRC} Info.h)

# When Kokkos targets the CUDA backend it compiles via nvcc. Allow calling
# constexpr __host__ functions (e.g. std::array::operator[]) from device code.
if (Kokkos_ENABLE_CUDA)
    target_compile_options(${PROJECT_NAME}_lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>)
endif ()
