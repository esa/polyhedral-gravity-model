include(FetchContent)

message(STATUS "Setting up spdlog")
set(SPDLOG_VERSION 1.14.1)

# Known Issue:
# If you install spdlog@1.14.1 via homebrew on ARM macOS, CMake will find spdlog
# However, there is a version mismatch between the `fmt` library installed as dependency, and the one actually
# being required leading to a linking error (i.e. missing symbols) while compiling!
find_package(spdlog ${SPDLOG_VERSION} QUIET)

if(${spdlog_FOUND})
    message(STATUS "Found existing spdlog Library: ${spdlog_DIR}")
else()
    message(STATUS "Using Spdlog Library from GitHub Release ${SPDLOG_VERSION}")
    FetchContent_Declare(spdlog
            GIT_REPOSITORY https://github.com/gabime/spdlog.git
            GIT_TAG v${SPDLOG_VERSION}
    )
    set(SPDLOG_BUILD_EXAMPLE OFF CACHE BOOL "" FORCE)
    set(SPDLOG_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(SPDLOG_INSTALL OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(spdlog)
endif()