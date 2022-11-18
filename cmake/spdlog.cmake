include(FetchContent)

message(STATUS "Setting up spdlog")

find_package(spdlog 1.11.0 QUIET)

if (${spdlog_FOUND})

    message(STATUS "Using existing spdlog installation")

else()

    message(STATUS "Using spdlog from git repository")

    #Fetches the version 1.11.0 for spdlog
    FetchContent_Declare(spdlog
            GIT_REPOSITORY https://github.com/gabime/spdlog.git
            GIT_TAG v1.11.0
            )

    # Disable stuff we don't need
    option(SPDLOG_BUILD_EXAMPLE "" OFF)
    option(SPDLOG_BUILD_TESTS "" OFF)
    option(SPDLOG_INSTALL "" OFF)

    FetchContent_MakeAvailable(spdlog)

    # Disable warnings from the library target
    target_compile_options(spdlog PRIVATE -w)
    # Disable warnings from included headers
    get_target_property(propval spdlog INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(spdlog SYSTEM PUBLIC "${propval}")

endif()