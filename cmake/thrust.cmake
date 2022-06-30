include(FetchContent)

message(STATUS "Setting up thrust")

find_package(Thrust 1.16.0 QUIET)

if (${Thrust_FOUND})

    message(STATUS "Using local thrust installation")

else()
    message(STATUS "Using thrust from git repository")
    # Fetches the version 1.16.0 of the official NVIDIA Thrust repository
    FetchContent_Declare(thrust
            GIT_REPOSITORY https://github.com/NVIDIA/thrust.git
            GIT_TAG 1.16.0
            )

    FetchContent_GetProperties(thrust)
    if(NOT thrust_POPULATED)
        # Fetch the content using previously declared details
        FetchContent_Populate(thrust)

        # Set custom variables, policies, etc.
        # Disable stuff not needed
        set(THRUST_ENABLE_HEADER_TESTING "OFF")
        set(THRUST_ENABLE_TESTING "OFF")
        set(THRUST_ENABLE_EXAMPLES "OFF")

        # Set standard CPP Dialect to 17 (default of thrust would be 14)
        set(THRUST_CPP_DIALECT 17)


        if(EXISTS ${thrust_SOURCE_DIR}/thrust/cmake/thrust-config-version.cmake)
            message(STATUS "Found problem causing thrust file (1)")
            file(READ ${thrust_SOURCE_DIR}/thrust/cmake/thrust-config-version.cmake THRUST_CONFIG_VERSION_FILE)
            string(REPLACE
                    "file(READ \"\${_THRUST_VERSION_INCLUDE_DIR}/thrust/version.h\" THRUST_VERSION_HEADER)"
                    "file(READ \"${thrust_SOURCE_DIR}/thrust/version.h\" THRUST_VERSION_HEADER)"
                    THRUST_CONFIG_VERSION_FILE "${THRUST_CONFIG_VERSION_FILE}")

            file(WRITE ${thrust_SOURCE_DIR}/thrust/cmake/thrust-config-version.cmake
                    "${THRUST_CONFIG_VERSION_FILE}"
                    )
            message(STATUS "Modified the the problem causing file (1)")
        else()
            message(STATUS "Problem causing thrust file (1) not found!")
        endif()

        if(EXISTS ${thrust_SOURCE_DIR}/thrust/cmake/thrust-config.cmake)
            message(STATUS "Found problem causing thrust file (2)")
            file(READ ${thrust_SOURCE_DIR}/thrust/cmake/thrust-config.cmake THRUST_CONFIG_FILE)
            string(REPLACE
                    "target_include_directories(_Thrust_Thrust INTERFACE \"\${_THRUST_INCLUDE_DIR}\")"
                    "target_include_directories(_Thrust_Thrust INTERFACE \"${thrust_SOURCE_DIR}\")"
                    THRUST_CONFIG_FILE "${THRUST_CONFIG_FILE}")

            file(WRITE ${thrust_SOURCE_DIR}/thrust/cmake/thrust-config.cmake
                    "${THRUST_CONFIG_FILE}"
                    )
            message(STATUS "Modified the the problem causing file (2)")
        else()
            message(STATUS "Problem causing thrust file (2) not found!")
        endif()

        # Bring the populated content into the build
        add_subdirectory(${thrust_SOURCE_DIR} ${thrust_BINARY_DIR})

    endif()
endif()
