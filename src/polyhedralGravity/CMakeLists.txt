# Including the actual Polyhedral source files
file(GLOB_RECURSE SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# The OpenCL backend is only compiled in when it is enabled
if (NOT POLYHEDRAL_GRAVITY_ENABLE_OPENCL)
    list(FILTER SRC EXCLUDE REGEX "/opencl/")
endif ()

add_library(${PROJECT_NAME}_lib OBJECT ${SRC} Info.h)

if (POLYHEDRAL_GRAVITY_ENABLE_OPENCL)
    # Plain signature, matching how src/CMakeLists.txt links this target
    target_link_libraries(${PROJECT_NAME}_lib polyhedral_gravity_opencl)
    # The generated kernel headers have to exist before the OpenCL sources are compiled
    add_dependencies(${PROJECT_NAME}_lib polyhedral_gravity_opencl_kernels)
endif ()
