function(polyhedral_gravity_parse_version OUTPUT_VAR)
    # Read the content of the given header file
    file(READ "${PROJECT_SOURCE_DIR}/src/polyhedralGravity/Version.h" HEADER_CONTENTS)
    # Extract the version using regex
    string(REGEX MATCH "constexpr std::string_view POLYHEDRAL_GRAVITY_VERSION = \"([0-9]+\\.[0-9]+\\.[0-9]+[a-zA-Z0-9]*)\"" VERSION_MATCH "${HEADER_CONTENTS}")
    # Set the output variable to the matched version group
    if(CMAKE_MATCH_1)
        set(${OUTPUT_VAR} "${CMAKE_MATCH_1}" PARENT_SCOPE)
    else()
        message(FATAL_ERROR "Failed to parse POLYHEDRAL_GRAVITY_VERSION from '${HEADER_FILE}'")
    endif()
endfunction()