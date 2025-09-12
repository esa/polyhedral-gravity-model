# Get the Git information
get_git_commit_hash(POLYHEDRAL_GRAVITY_COMMIT_HASH)
is_git_working_tree_clean(POLYHEDRAL_GRAVITY_WORKING_TREE)
# The Regex of the setup.py captures everything inside the quotes!
set(POLYHEDRAL_GRAVITY_VERSION "3.3.1rc1")
set(${PROJECT_VERSION} ${POLYHEDRAL_GRAVITY_VERSION})

# Append "-modified" to the commit hash if the working tree is not clean
if (NOT ${POLYHEDRAL_GRAVITY_WORKING_TREE})
    set(POLYHEDRAL_GRAVITY_COMMIT_HASH "${POLYHEDRAL_GRAVITY_COMMIT_HASH}+modified")
endif ()

# Configure the output header file
configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/src/polyhedralGravity/Info.h.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/polyhedralGravity/Info.h"
)