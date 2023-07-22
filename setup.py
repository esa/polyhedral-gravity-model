import os
import re
import subprocess
import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

# ---------------------------------------------------------------------------------
# Modify these variables to customize the build of the polyhedral gravity interface
# All variables can be overwritten by setting equally named env variables
# ---------------------------------------------------------------------------------
# The generator used for CMake
CMAKE_GENERATOR = "Ninja"
# The other CMake options
CMAKE_OPTIONS = {
    # The Build Type (Should be release!)
    "CMAKE_BUILD_TYPE": "Release",
    # Modify to change the parallelization (Default value: TBB)
    "POLYHEDRAL_GRAVITY_PARALLELIZATION": "TBB",
    # Default value (INFO=2)
    "LOGGING_LEVEL": 2,
    # Default value (OFF)
    "USE_LOCAL_TBB": "OFF",
    # Not required for the python interface (--> OFF)
    "BUILD_POLYHEDRAL_GRAVITY_DOCS": "OFF",
    # Not required for the python interface (--> OFF)
    "BUILD_POLYHEDRAL_GRAVITY_TESTS": "OFF",
    # Not required for the python interface (--> OFF)
    "BUILD_POLYHEDRAL_GRAVITY_LIBRARY": "OFF",
    # Not required for the python interface (--> OFF)
    "BUILD_POLYHEDRAL_GRAVITY_EXECUTABLE": "OFF",
    # Should be of course ON!
    "BUILD_POLYHEDRAL_GRAVITY_PYTHON_INTERFACE": "ON",
    # Build static libs by default (On conda-forge we build shared libs by setting this to ON)
    "BUILD_SHARED_LIBS": "OFF"
}
# ---------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
# The following is adapted from https://github.com/pybind/cmake_example/blob/master/setup.py
# -----------------------------------------------------------------------------------------

# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection & inclusion of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
        ]
        build_args = []

        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        # (no specified CMAKE_ARGS not covered by CMAKE_OPTIONS)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        # Get the path/to/conda/env and include it in the cmake prefix path
        prefix = os.environ.get("PREFIX", "")
        if prefix:
            cmake_args += [
                f"-DCMAKE_PREFIX_PATH={prefix}"
            ]

        # Use the values specified in the script if no env variable (like on conda) overwrites it
        for option_name, option_default_value in CMAKE_OPTIONS.items():
            final_value = os.environ.get(option_name, option_default_value)
            cmake_args += [
                f"-D{option_name}={final_value}"
            ]

        # Sets the CMake Generator if specified (this is separate from the other variables since it is given to
        # CMake vie the -G prefix
        final_generator = os.environ.get("CMAKE_GENERATOR", CMAKE_GENERATOR)
        cmake_args += [
            f"-G{final_generator}"
        ]

        # MSVC special cases
        if self.compiler.compiler_type == "msvc":
            # Single config generators are handled "normally"
            single_config = any(x in CMAKE_GENERATOR for x in {"NMake", "Ninja"})
            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{CMAKE_OPTIONS['CMAKE_BUILD_TYPE'].upper()}={extdir}"
                ]
                build_args += ["--config", CMAKE_OPTIONS['CMAKE_BUILD_TYPE']]

        # MacOS special cases
        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += [f"-j{self.parallel}"]

        build_temp = os.path.join(self.build_temp, ext.name)
        if not os.path.exists(build_temp):
            os.makedirs(build_temp)

        # Enables log messages for thrust
        cmake_args += ["--log-level=VERBOSE"]

        # Call CMake and build the project
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=build_temp)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=build_temp)


# -----------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# Modify these entries to customize the package metadata of the polyhedral gravity
# --------------------------------------------------------------------------------
setup(
    name="polyhedral_gravity",
    version="2.0",
    author="Jonas Schuhmacher",
    author_email="jonas.schuhmacher@tum.de",
    description="Package to compute full gravity tensor of a given constant density polyhedron for arbitrary points",
    long_description="""
        The package polyhedral_gravity provides a simple to use interface for the evaluation of the full gravity
        tensor of a constant density polyhedron at given computation points. It is based on a fast, parallelized
        backbone in C++ capable of evaluating the gravity at thousands of computation points in the fraction of a second.
    """,
    ext_modules=[CMakeExtension("polyhedral_gravity")],
    cmdclass={"build_ext": CMakeBuild},
    license_file="LICENSE",
    zip_safe=False,
    python_requires=">=3.6",
    include_package_data=True,
)
# --------------------------------------------------------------------------------
