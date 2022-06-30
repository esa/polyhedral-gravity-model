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

        # Original Content of the problematic file
        ## Parse version information from version.h:
        ##unset(_THRUST_VERSION_INCLUDE_DIR CACHE) # Clear old result to force search
        ##find_path(_THRUST_VERSION_INCLUDE_DIR thrust/version.h
        ##  NO_DEFAULT_PATH # Only search explicit paths below:
        ##  PATHS
        ##    "${CMAKE_CURRENT_LIST_DIR}/../.."            # Source tree
        ##)
        ##set_property(CACHE _THRUST_VERSION_INCLUDE_DIR PROPERTY TYPE INTERNAL)

        if(EXISTS ${thrust_SOURCE_DIR}/thrust/cmake/thrust-header-search.cmake)
            message(STATUS "Found problem causing thrust file!")

            file(WRITE ${thrust_SOURCE_DIR}/thrust/cmake/thrust-header-search.cmake
                    "
# Parse version information from version.h:
unset(_THRUST_VERSION_INCLUDE_DIR CACHE) # Clear old result to force search
find_path(_THRUST_VERSION_INCLUDE_DIR thrust/version.h
  # NO_DEFAULT_PATH # Only search explicit paths below: COMMENTED BY POLYHEDRAL_GRAVITY_CMAKE
  PATHS
    \"\${CMAKE_CURRENT_LIST_DIR}/../..\"            # Source tree
)
set_property(CACHE _THRUST_VERSION_INCLUDE_DIR PROPERTY TYPE INTERNAL)
                    ")
            message(STATUS "Modified the the problem causing file!")
        else()
            message(STATUS "Problem causing thrust file not found!")
        endif()

        # Bring the populated content into the build
        add_subdirectory(${thrust_SOURCE_DIR} ${thrust_BINARY_DIR})

    endif()
endif()
