# polyhedral-gravity-model

![Build and Test](https://github.com/schuhmaj/polyhedral-gravity-model-cpp/actions/workflows/ctest.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/polyhedral-gravity-model-cpp/badge/?version=stable)](https://polyhedral-gravity-model-cpp.readthedocs.io/en/stable/?badge=stable)
![GitHub](https://img.shields.io/github/license/esa/polyhedral-gravity-model)

![PyPI](https://img.shields.io/pypi/v/polyhedral-gravity)
![Static Badge](https://img.shields.io/badge/platform-linux--64_%7C_win--64_%7C_osx--64_%7C_linux--arm64_%7C_osx--arm64-lightgrey)
![PyPI - Downloads](https://img.shields.io/pypi/dm/polyhedral-gravity)

![Conda](https://img.shields.io/conda/v/conda-forge/polyhedral-gravity-model)
![Conda](https://img.shields.io/conda/pn/conda-forge/polyhedral-gravity-model)
![Conda](https://img.shields.io/conda/dn/conda-forge/polyhedral-gravity-model)

## Table of Contents

- [References](#references)
- [Documentation & Examples](#documentation--examples)
  - [Overview](#overview)
  - [Minimal Python Example](#minimal-python-example)
  - [Minimal C++ Example](#minimal-c-example)
- [Installation](#installation)
  - [With conda](#with-conda)
  - [With pip](#with-pip)
  - [From source](#from-source)
- [Miscellaneous](#miscellaneous)
  - [Building the C++ Library & Executable](#building-the-c-library--executable)
  - [Supported Polyhedron Source Files (Python/ C++)](#supported-polyhedron-source-files-python-c)
  - [The C++ Executable](#the-c-executable)
    - [Config File](#config-file)
    - [Output](#output)
- [Testing](#testing)
- [Contributing](#contributing)

## References

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

## Documentation & Examples

### Overview

Some exemplary results and plots are stored in the
[jupyter notebook](script/polyhedral-gravity.ipynb).
It also provides a good introduction to the application of
the python interface.

The full extensive documentation can be found
on [readthedocs](https://polyhedral-gravity-model-cpp.readthedocs.io/en/stable/).

### Minimal Python Example

The following example shows how to use the python interface to compute the gravity
around a cube:

```python
import numpy as np
import polyhedral_gravity

# We define the cube as a polyhedron with 8 vertices and 12 triangular faces
# The density is set to 1.0
cube_vertices = np.array(
    [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]
)
cube_faces = np.array(
    [[1, 3, 2], [0, 3, 1], [0, 1, 5], [0, 5, 4], [0, 7, 3], [0, 4, 7],
    [1, 2, 6], [1, 6, 5], [2, 3, 6], [3, 7, 6], [4, 5, 6], [4, 6, 7]]
)
cube_density = 1.0
computation_point = np.array([[0, 0, 0]])
```

The simplest way to compute the gravity is to use the `evaluate` function:

```python
potential, acceleration, tensor = polyhedral_gravity.evaluate(
    polyhedral_source=(cube_vertices, cube_faces),
    density=cube_density,
    computation_points=computation_point
    parallel=True
)
```

The more advanced way is to use the `GravityEvaluable` class. It caches the
internal data strcuture and properties which can be reused for multiple
evaluations. This is especially useful if you want to compute the gravity
for multiple computation points, but don't know the "future points" in advance.

```python
evaluable = model.GravityEvaluable(
    polyhedral_source=(cube_vertices, cube_faces),
    density=cube_density
)
potential, acceleration, tensor = evaluable(computation_point, parallel=True)
```

In case you want to hand over the polyhedron via a supported file format,
just replace the `polyhedral_source` argument with *a list of strings*,
where each string is the path to a supported file format.

### Minimal C++ Example

The following example shows how to use the C++ library to compute the gravity.
It works analogously to the Python example above.

```cpp
// Defining the input like above in the Python example
std::vector<std::array<double, 3>> vertices = ...
std::vector<std::array<size_t, 3>> faces = ...
// The polyhedron is defined by its vertices and faces
Polyhedron polyhedron{vertices, faces};
double density = 1.0;
std::vector<std::array<double, 3>> points = ...
std::array<double, 3> point = points[0];
bool parallel = true;
```

The C++ library provides also two ways to compute the gravity. Via
the free function `evaluate`...

```cpp
const auto[pot, acc, tensor] = GravityModel::evaluate(polyhedron, density, point);
```

... or via the `GravityEvaluable` class.

```cpp
// Instantiation of the GravityEvaluable object
GravityEvaluable evaluable{polyhedron, density};

// From now, we can evaluate the gravity model for any point with
const auto[potential, acceleration, tensor] = evaluable(point, parallel);
// or for multiple points with
const auto results = evaluable(points);
```

## Installation

### With conda

The python interface can be easily installed with
[conda](https://anaconda.org/conda-forge/polyhedral-gravity-model):

```bash
conda install -c conda-forge polyhedral-gravity-model
```

### With pip

As a second option, you can also install the python interface with pip.

```bash
pip install polyhedral-gravity
```

Binaries for the most common platforms are available on PyPI including
Windows, Linux and macOS. For macOS and Linux, binaries for
`x86_64` and `aarch64` are provided.
In case `pip` uses the source distribution, please make sure that
you have a C++17 capable compiler, CMake and ninja-build installed.

### From source

The project uses the following dependencies,
all of them are **automatically** set-up via CMake:

- GoogleTest (1.13.0 or compatible), only required for testing
- spdlog (1.11.0 or compatible), required for logging
- tetgen (1.6 or compatible), required for I/O
- yaml-cpp (0.7.0 or compatible), required for I/O
- thrust (2.1.0 or compatible), required for parallelization and utility
- xsimd (11.1.0 or compatible), required for vectorization of the `atan(..)`
- pybind11 (2.10.4 or compatible), required for the Python interface, but not the C++ standalone

The module will be build using a C++17 capable compiler,
CMake and ninja-build. Just execute the following command in
the repository root folder:

```bash
pip install .
```

To modify the build options (like parallelization) have a look
at the `setupy.py` and the [next paragraph](#building-the-c-library--executable). The options
are modified by setting the environment variables before executing
the `pip install .` command, e.g.:

```bash
export POLYHEDRAL_GRAVITY_PARALLELIZATION="TBB"
pip install .
```

(Optional: For a faster build you can install all dependencies available
for your system in your local python environment. That way, they
won't be fetched from GitHub.)

## Miscellaneous

### Building the C++ Library & Executable

The program is build by using CMake. So first make sure that you installed
CMake and then follow these steps:

```bash
mkdir build
cd build
cmake .. <options>
cmake --build .
```

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

```bash
cmake .. -POLYHEDRAL_GRAVITY_PARALLELIZATION="TBB"
```

### Supported Polyhedron Source Files (Python/ C++)

The implementation supports multiple common mesh formats for
the polyhedral source. These fromats are supported by the C++ library
and the Python interface.
These include:

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

### The C++ Executable

After the build, the gravity model can be run by executing:

```bash
    ./polyhedralGravity <YAML-Configuration-File>
```

where the YAML-Configuration-File contains the required parameters.
Examples for Configuration Files and Polyhedral Source Files can be
found in this repository in the folder `/example-config/`.

#### Config File

The configuration should look similar to the given example below.
It is required to specify the source-files of the polyhedron's mesh (more info
about the supported file in the [previous paragraph](#supported-polyhedron-source-files-python-c)), the density
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

#### Output

The calculation outputs the following parameters for every Computation Point _P_:

|             Name             |      Unit       |                              Comment                              |
|:----------------------------:|:---------------:|:-----------------------------------------------------------------:|
|              V               | m^2/s^2 or J/kg |           The potential or also called specific energy            |
|          Vx, Vy, Vz          |      m/s^2      | The gravitational accerleration in the three cartesian directions |
| Vxx, Vyy, Vzz, Vxy, Vxz, Vyz |      1/s^2      |   The spatial rate of change of the gravitational accleration    |

## Testing

The project uses GoogleTest for testing. In oder to execute those
tests just execute the following command in the build directory:

```bash
ctest
```

For the Python test suite, please execute the following command in the repository root folder:

```bash
pytest
```

## Contributing

We are happy to accept contributions to the project in the form of
suggestions, bug reports and pull requests. Please have a look at
the [contributing guidelines](CONTRIBUTING.md) for more information.
