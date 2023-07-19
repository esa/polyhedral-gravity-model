# polyhedral-gravity-model

![Build and Test](https://github.com/schuhmaj/polyhedral-gravity-model-cpp/actions/workflows/ctest.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/polyhedral-gravity-model-cpp/badge/?version=latest)](https://polyhedral-gravity-model-cpp.readthedocs.io/en/latest/?badge=latest)

This code is a validated implementation in C++17 of the Polyhedral Gravity Model
by Tsoulis et al.. It was created in a collaborative project between
TU Munich and ESA's Advanced Concepts Team. Please refer to the
[project report](https://mediatum.ub.tum.de/1695208)
for extensive information about the theoretical background, related work,
implementation & design decisions, application, verification,
and runtime measurements of the presented code.


The implementation is based on the
paper [Tsoulis, D., 2012. Analytical computation of the full gravity tensor of a homogeneous arbitrarily shaped polyhedral source using line integrals. Geophysics, 77(2), pp.F1-F11.](http://dx.doi.org/10.1190/geo2010-0334.1)
and its corresponding [implementation in FORTRAN](https://software.seg.org/2012/0001/index.html).

Supplementary details can be found in the more recent
paper [TSOULIS, Dimitrios; GAVRIILIDOU, Georgia. A computational review of the line integral analytical formulation of the polyhedral gravity signal. Geophysical Prospecting, 2021, 69. Jg., Nr. 8-9, S. 1745-1760.](https://doi.org/10.1111/1365-2478.13134)
and its corresponding [implementation in MATLAB](https://github.com/Gavriilidou/GPolyhedron),
which is strongly based on the former implementation in FORTRAN.

### Results and Plots

Some exemplary results and plots are stored in the
[jupyter notebook](script/polyhedral-gravity.ipynb).
It also provides a good introduction to the application of
the python interface.

## Documentation (readthedocs)

The full extensive documentation can be found
on [readthedocs](https://polyhedral-gravity-model-cpp.readthedocs.io/en/latest/).

## Requirements

The project uses the following dependencies,
all of them are **automatically** set-up via CMake:

- GoogleTest (1.13.0 or compatible), only required for testing
- spdlog (1.11.0 or compatible), required for logging
- tetgen (1.6 or compatible), required for I/O
- yaml-cpp (0.7.0 or compatible), required for I/O
- thrust (2.1.0 or compatible), required for parallelization and utility
- xsimd (11.1.0 or compatible), required for vectorization of the `atan(..)`
- pybind11 (2.10.4 or compatible), required for the Python interface, but not the C++ standalone

## Python interface

### conda

The python interface can be easily installed with
[conda](https://anaconda.org/conda-forge/polyhedral-gravity-model):

    conda install -c conda-forge polyhedral-gravity-model

This is currently only supported for `x86-64` systems since
one of the dependencies is not available on conda for `aarch64`.
However, building from source with `pip` can also be done
on `aarch64` as shown below.

### pip

Use pip to install the python interface in your local python runtime.
The module will be build using CMake and the using the above
requirements. Just execute in repository root:

    pip install .

To modify the build options (like parallelization) have a look
at the `setupy.py` and the [next paragraph](#build-c).
(Optional: For a faster build you can install all dependencies available
for your system in your local python environment. That way, they
won't be fetched from GitHub.)

## Build C++

The program is build by using CMake. So first make sure that you installed
CMake and then follow these steps:

    mkdir build
    cd build
    cmake .. <options>
    cmake --build .

The following options are available:

|                             Name (Default) | Options                                                                                    |
|-------------------------------------------:|:-------------------------------------------------------------------------------------------|
| POLYHEDRAL_GRAVITY_PARALLELIZATION (`CPP`) | `CPP` = Serial Execution / `OMP` or `TBB` = Parallel Execution with OpenMP or Intel\'s TBB |
|                        LOGGING_LEVEL (`2`) | `0` = TRACE/ `1` = DEBUG/ `2` = INFO / `3` = WARN/ `4` = ERROR/ `5` = CRITICAL/ `6` = OFF  |
|                      USE_LOCAL_TBB (`OFF`) | Use a local installation of `TBB` instead of setting it up via `CMake`                     |
|      BUILD_POLYHEDRAL_GRAVITY_DOCS (`OFF`) | Build this documentation                                                                   |
|      BUILD_POLYHEDRAL_GRAVITY_TESTS (`ON`) | Build the Tests                                                                            |
|   BUILD_POLYHEDRAL_PYTHON_INTERFACE (`ON`) | Build the Python interface                                                                 |

During testing POLYHEDRAL_GRAVITY_PARALLELIZATION=`TBB` has been the most performant.
It is further not recommend to change the LOGGING_LEVEL to something else than `INFO=2`.

The recommended CMake command would look like this (we only need to change `PARALLELIZATION_DEVICE`, since
the defaults of the others are already correctly set):

    cmake .. -POLYHEDRAL_GRAVITY_PARALLELIZATION="TBB"

## Execution C++

### Overview

After the build, the gravity model can be run by executing:

    ./polyhedralGravity <YAML-Configuration-File>

where the YAML-Configuration-File contains the required parameters.
Examples for Configuration Files and Polyhedral Source Files can be
found in this repository in the folder `/example-config/`.

### Config File

The configuration should look similar to the given example below.
It is required to specify the source-files of the polyhedron's mesh (more info
about the supported file in the [next paragraph](#polyhedron-source-files)), the density
of the polyhedron, and the wished computation points where the
gravity tensor shall be computed.
Further one must specify the name of the .csv output file.

````yaml
---
gravityModel:
  input:
    polyhedron: #polyhedron source-file(s)
      - "../example-config/data/tsoulis.node"   # .node contains the vertices
      - "../example-config/data/tsoulis.face"   # .face contains the triangular faces
    density: 2670.0                             # constant density in [kg/m^3]
    points: # Location of the computation point(s) P
      - [ 0, 0, 0 ]                             # Here it is situated at the origin
    check_mesh: true                            # Fully optional, enables input checking (not given: false)
  output:
    filename: "gravity_result.csv"              # The name of the output file 

````

### Polyhedron Source Files

The implementation supports multiple common mesh formats for
the polyhedral source. These include:

|     File Suffix     |                        Name                        | Comment                                                                                                                                          |
|:-------------------:|:--------------------------------------------------:|--------------------------------------------------------------------------------------------------------------------------------------------------|
| `.node` and `.face` |                   TetGen's files                   | These two files need to be given as a pair to the input. [Documentation of TetGen's files](https://wias-berlin.de/software/tetgen/fformats.html) |
|       `.mesh`       |                 Medit's mesh files                 | Single file containing every needed mesh information.                                                                                            |
|       `.ply`        | The Polygon File format/ Stanfoard Triangle format | Single file containing every needed mesh information. Blender File Format.                                                                       |
|       `.off`        |                 Object File Format                 | Single file containing every needed mesh information.                                                                                            |
|       `.stl`        |              Stereolithography format              | Single file containing every needed mesh information. Blender File Format.                                                                       |                                         

**Notice!** Only the ASCII versions of those respective files are supported! This is especially
important for e.g. the `.ply` files which also can be in binary format.

Good tools to convert your Polyhedron to a supported format (also for interchanging
ASCII and binary format) are e.g.:

- [Meshio](https://github.com/nschloe/meshio) for Python
- [OpenMesh](https://openmesh-python.readthedocs.io/en/latest/readwrite.html) for Python

The vertices in the input mesh file must be ordered so that the plane unit normals point outwards of the polyhedron for
every face.
One can use the program input-checking procedure to ensure the correct format. This method is activated via the
corresponding configuration option and uses the Möller–Trumbore intersection algorithm. Notice that this algorithm is a
quadratic complexity, so the check should only be utilized in case of uncertainty.

### Output

The calculation outputs the following parameters for every Computation Point _P_:

|             Name             |      Unit       |                              Comment                              |
|:----------------------------:|:---------------:|:-----------------------------------------------------------------:|
|              V               | m^2/s^2 or J/kg |           The potential or also called specific energy            |
|          Vx, Vy, Vz          |      m/s^2      | The gravitational accerleration in the three cartesian directions |
| Vxx, Vyy, Vzz, Vxy, Vxz, Vyz |      1/s^2      |   The spatial rate of change of the gravitational accleration    |

## Testing C++

The project uses GoogleTest for testing. In oder to execute those
tests just execute the following command in the build directory:

    ctest
