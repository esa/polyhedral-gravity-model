include(FetchContent)

message(STATUS "Setting up tetgen")

# IMPORTANT NOTE
# We do ot use findPackage here, as we modify the one source file slightly to suppress output via stdout!!!

#Fetches the version 1.6 for tetgen
FetchContent_Declare(tetgen
        GIT_REPOSITORY https://github.com/libigl/tetgen.git
        GIT_TAG 4f3bfba3997f20aa1f96cfaff604313a8c2c85b6 # release 1.6
        )

FetchContent_GetProperties(tetgen)
if(NOT tetgen_POPULATED)
    FetchContent_Populate(tetgen)

    # This is ugly, but functional
    # We modify one source file in order to prevent console output from the printf function
    # without polluting the global include path (which would cause to more uglier issues)
    if(NOT EXISTS ${tetgen_SOURCE_DIR}/tetgen_mod.cxx)
        message(STATUS "Creating modified tetgen.cxx in order to prevent console output from library")
        file(READ ${tetgen_SOURCE_DIR}/tetgen.cxx TETGEN_CXX)

        string(REPLACE
                "#include \"tetgen.h\""
                "#include \"tetgen.h\"\n#define printf(fmt, ...) (0)\n"
                TETGEN_CXX "${TETGEN_CXX}")

        file(WRITE ${tetgen_SOURCE_DIR}/tetgen_mod.cxx
                "${TETGEN_CXX}"
                )
    else()
        message(STATUS "A modified tetgen.cxx already exists! It is assumed that is the correct one disabling output")
    endif()

    add_library(tetgen STATIC
            ${tetgen_SOURCE_DIR}/tetgen_mod.cxx
            ${tetgen_SOURCE_DIR}/predicates.cxx
            )

    target_compile_definitions(tetgen PRIVATE -DTETLIBRARY)

    target_include_directories(tetgen INTERFACE "${tetgen_SOURCE_DIR}")

endif()

# Disable warnings from the library target
target_compile_options(tetgen PRIVATE -w)