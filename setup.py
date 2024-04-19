import os
import re
import subprocess
import shutil
import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

# ---------------------------------------------------------------------------------
# Modify these variables to customize the build of the polyhedral gravity interface
# All variables can be overwritten by setting equally named env variables
# ---------------------------------------------------------------------------------
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


def is_ninja_installed():
    """Returns true if ninja build tool is installed. False otherwise."""
    return shutil.which("ninja") is not None


def get_cmake_generator():
    """Returns the CMake generator if specified as environment variable.
    If not, returns "Ninja" if ninja build is installed.
    Otherwise, none is returned.
    """
    os_env_generator = os.environ.get("CMAKE_GENERATOR", None)
    if os_env_generator is not None:
        return os_env_generator
    elif is_ninja_installed():
        return "Ninja"
    else:
        return None

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

        # Disable availability of standard libc++ on macOS if requested
        if os.environ.get("_LIBCPP_DISABLE_AVAILABILITY"):
            cmake_args += [
                "-D_LIBCPP_DISABLE_AVAILABILITY=ON"
            ]

        # Sets the CMake Generator if specified (this is separate from the other variables since it is given to
        # CMake vie the -G prefix
        cmake_generator = get_cmake_generator()
        if cmake_generator is not None:
            cmake_args += [
                f"-G{cmake_generator}"
            ]

        # MSVC special cases
        if self.compiler.compiler_type == "msvc":
            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})
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


picture_in_readme = '''<p align="center">
  <img src="paper/figures/eros_010.png" width="50%">
  <br>
  <em>
    <a href="https://github.com/gomezzz/geodesyNets/blob/master/3dmeshes/eros_lp.pk">Mesh of (433) Eros</a> with 739 vertices and 1474 faces
  </em>
</p>'''


# --------------------------------------------------------------------------------
# Modify these entries to customize the package metadata of the polyhedral gravity
# --------------------------------------------------------------------------------
setup(
    name="polyhedral_gravity",
    version="3.0rc4",
    author="Jonas Schuhmacher",
    author_email="jonas.schuhmacher@tum.de",
    description="Package to compute full gravity tensor of a given constant density polyhedron for arbitrary points "
                "according to the geodetic convention",
    long_description=open("README.md").read().replace(picture_in_readme, ''),
    long_description_content_type="text/markdown",
    ext_modules=[CMakeExtension("polyhedral_gravity")],
    cmdclass={"build_ext": CMakeBuild},
    license="GPLv3",
    license_file="LICENSE",
    zip_safe=False,
    python_requires=">=3.6",
    include_package_data=True,
    project_urls={
        "Homepage": "https://github.com/esa/polyhedral-gravity-model",
        "Source": "https://github.com/esa/polyhedral-gravity-model",
        "Documentation": "https://esa.github.io/polyhedral-gravity-model/",
        "Issues": "https://github.com/esa/polyhedral-gravity-model/issues",
        "Changelog": "https://github.com/esa/polyhedral-gravity-model/releases",
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
# --------------------------------------------------------------------------------
