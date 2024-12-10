function(get_git_commit_hash OUTPUT_VAR)
    # Ensure Git is available
    find_package(Git QUIET REQUIRED)

    # Run a Git command to get the current commit hash
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT_HASH
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_VARIABLE GIT_ERROR
            ERROR_STRIP_TRAILING_WHITESPACE
    )

    # Check if the Git command was successful
    if (NOT GIT_COMMIT_HASH OR GIT_ERROR)
        message(WARNING "Failed to retrieve Git commit hash: ${GIT_ERROR}")
        set(${OUTPUT_VAR} "UNKNOWN" PARENT_SCOPE)
    else()
        # Pass the hash back to the calling scope
        set(${OUTPUT_VAR} "${GIT_COMMIT_HASH}" PARENT_SCOPE)
    endif()
endfunction()