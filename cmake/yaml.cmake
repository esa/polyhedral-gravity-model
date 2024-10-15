include(FetchContent)

message(STATUS "Setting up yaml-cpp")
set(YAML_CPP_VERSION 0.8.0)

find_package(yaml-cpp ${YAML_CPP_VERSION} QUIET)

if (${yaml-cpp_FOUND})
    message(STATUS "Found existing yaml-cpp library: ${yaml-cpp_DIR}")
else()
    message(STATUS "Using yaml-cpp from GitHub Release ${YAML_CPP_VERSION}")
    FetchContent_Declare(yaml-cpp
            GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
            GIT_TAG ${YAML_CPP_VERSION}
            )
    # Disable everything we don't need
    set(YAML_CPP_BUILD_TESTS OFF CACHE INTERNAL "")
    set(YAML_CPP_BUILD_CONTRIB OFF CACHE INTERNAL "")
    set(YAML_CPP_BUILD_TOOLS OFF CACHE INTERNAL "")
    set(YAML_CPP_FORMAT_SOURCE OFF CACHE INTERNAL "")
    FetchContent_MakeAvailable(yaml-cpp)
    # Disable warnings from the library target
    target_compile_options(yaml-cpp PRIVATE -w)
    # Disable warnings from included headers
    get_target_property(propval yaml-cpp INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(yaml-cpp SYSTEM PUBLIC "${propval}")
endif()
