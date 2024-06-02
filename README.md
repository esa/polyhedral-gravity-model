# polyhedral-gravity-model

[![DOI](https://joss.theoj.org/papers/10.21105/joss.06384/status.svg)](https://doi.org/10.21105/joss.06384)
![GitHub](https://img.shields.io/github/license/esa/polyhedral-gravity-model)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/esa/polyhedral-gravity-model/.github%2Fworkflows%2Fbuild-and-test.yml?logo=GitHub%20Actions&label=Build%20and%20Test)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/esa/polyhedral-gravity-model/.github%2Fworkflows%2Fdocs.yml?logo=GitBook&label=Documentation)

![PyPI](https://img.shields.io/pypi/v/polyhedral-gravity)
![Static Badge](https://img.shields.io/badge/platform-linux--64_%7C_win--64_%7C_osx--64_%7C_linux--arm64_%7C_osx--arm64-lightgrey)
![PyPI - Downloads](https://img.shields.io/pypi/dm/polyhedral-gravity)

![Conda](https://img.shields.io/conda/v/conda-forge/polyhedral-gravity-model)
![Conda](https://img.shields.io/conda/pn/conda-forge/polyhedral-gravity-model)
![Conda](https://img.shields.io/conda/dn/conda-forge/polyhedral-gravity-model)

<p align="center">
  <img src="paper/figures/eros_010.png" width="50%">
  <br>
  <em>
    <a href="https://github.com/gomezzz/geodesyNets/blob/master/3dmeshes/eros_lp.pk">Mesh of (433) Eros</a> with 739 vertices and 1474 faces
  </em>
</p>

## Table of Contents

- [References](#references)
- [Documentation & Examples](#documentation--examples)
  - [Input & Output (C++ and Python)](#input--output-c-and-python)
  - [Minimal Python Example](#minimal-python-example)
  - [Minimal C++ Example](#minimal-c-example)
- [Installation](#installation)
  - [With conda](#with-conda)
  - [With pip](#with-pip)
  - [From source](#from-source)
- [C++ Library & Executable](#c-library--executable)
  - [Building the C++ Library & Executable](#building-the-c-library--executable)
  - [Running the C++ Executable](#running-the-c-executable)
- [Testing](#testing)
- [Contributing](#contributing)

## References

This code is a validated implementation in C++17 of the Polyhedral Gravity Model
by Tsoulis et al.. It was created in a collaborative project between
TU Munich and ESA's Advanced Concepts Team. Please refer to the
[project report](https://mediatum.ub.tum.de/doc/1695208/1695208.pdf)
for extensive information about the theoretical background, related work,
implementation & design decisions, application, verification,
and runtime measurements of the presented code.

The implementation is based on the
paper [Tsoulis, D., 2012. Analytical computation of the full gravity tensor of a homogeneous arbitrarily shaped polyhedral source using line integrals. Geophysics, 77(2), pp.F1-F11.](http://dx.doi.org/10.1190/geo2010-0334.1)
and its corresponding implementation in FORTRAN.

Supplementary details can be found in the more recent
paper [TSOULIS, Dimitrios; GAVRIILIDOU, Georgia. A computational review of the line integral analytical formulation of the polyhedral gravity signal. Geophysical Prospecting, 2021, 69. Jg., Nr. 8-9, S. 1745-1760.](https://doi.org/10.1111/1365-2478.13134)
and its corresponding [implementation in MATLAB](https://github.com/Gavriilidou/GPolyhedron),
which is strongly based on the former implementation in FORTRAN.

## Documentation & Examples

> [!NOTE]
> The [GitHub Pages](https://esa.github.io/polyhedral-gravity-model) of this project
contain the full extensive documentation of the C++ Library and Python Interface
as well as background on the gravity model and advanced settings not detailed here.

## Input & Output (C++ and Python)

### Input

The evaluation of the polyhedral gravity model requires the following parameters:

| Name                                                                       |
|----------------------------------------------------------------------------|
| Polyhedral Mesh (either as vertices & faces or as polyhedral source files) |
| Constant Density $\rho$                                                    |

The mesh and the constants density's unit must match.
Have a look the documentation to view the [supported mesh files](https://esa.github.io/polyhedral-gravity-model/quickstart/supported_input.html).

### Output

The calculation outputs the following parameters for every Computation Point *P*.
The units of the respective output depend on the units of the input parameters (mesh and density)!
Hence, if e.g. your mesh is in $km$, the density must match. Further, output units will be different accordingly.

|                            Name                            | Unit (if mesh in $[m]$ and $\rho$ in $[kg/m^3]$) |                              Comment                              |
|:----------------------------------------------------------:|:------------------------------------------------:|:-----------------------------------------------------------------:|
|                            $V$                             |       $\frac{m^2}{s^2}$ or $\frac{J}{kg}$        |           The potential or also called specific energy            |
|                    $V_x$, $V_y$, $V_z$                     |                 $\frac{m}{s^2}$                  | The gravitational accerleration in the three cartesian directions |
| $V_{xx}$, $V_{yy}$, $V_{zz}$, $V_{xy}$, $V_{xz}$, $V_{yz}$ |                 $\frac{1}{s^2}$                  |   The spatial rate of change of the gravitational accleration    |


>[!NOTE]
>This gravity model's output obeys to the geodesy and geophysics sign conventions.
Hence, the potential $V$ for a polyhedron with a mass $m > 0$ is defined as **positive**.
Accordingly, the accelerations are defined as $\textbf{g} = + \nabla V$.


### Minimal Python Example

The following example shows how to use the python interface to compute the gravity
around a cube:

```python
import numpy as np
from polyhedral_gravity import Polyhedron, GravityEvaluable, evaluate, PolyhedronIntegrity, NormalOrientation

# We define the cube as a polyhedron with 8 vertices and 12 triangular faces
# The polyhedron's normals point outwards (see below for checking this)
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
computation_point = np.array([0, 0, 0])
```

We first define a constant density Polyhedron from `vertices` and `faces`

```python
cube_polyhedron = Polyhedron(
  polyhedral_source=(cube_vertices, cube_faces),
  density=cube_density,
)
```

In case you want to hand over the polyhedron via a [supported file format](https://esa.github.io/polyhedral-gravity-model/quickstart/supported_input.html),
just replace the `polyhedral_source` argument with *a list of strings*,
where each string is the path to a supported file format, e.g. `polyhedral_source=["eros.node","eros.face"]` or `polyhedral_source=["eros.mesh"]`.

Continuing, the simplest way to compute the gravity is to use the `evaluate` function:

```python
potential, acceleration, tensor = evaluate(
  polyhedron=cube_polyhedron,
  computation_points=computation_point,
  parallel=True,
)
```

The more advanced way is to use the `GravityEvaluable` class. It caches the
internal data structure and properties which can be reused for multiple
evaluations. This is especially useful if you want to compute the gravity
for multiple computation points, but don't know the "future points" in advance.

```python
evaluable = GravityEvaluable(polyhedron=cube_polyhedron) # stores intermediate computation steps
potential, acceleration, tensor = evaluable(
  computation_points=computation_point,
  parallel=True,
)
# Any future evaluable call after this one will be faster
```

Note that the `computation_point` could also be (N, 3)-shaped array to compute multiple points at once.
In this case, the return value of `evaluate(..)` or an `GravityEvaluable` will
be a list of triplets comprising potential, acceleration, and tensor.

The gravity model requires that all the polyhedron's plane unit normals consistently
point outwards or inwards the polyhedron. You can specify this via the `normal_orientation`.
This property is - by default - checked when constructing the `Polyhedron`! So, don't worry, it
is impossible if not **explicitly** disabled to create an invalid `Polyhedron`.
You can disable/ enable this setting via the optional `integrity_check` flag and can even
automatically repair the ordering via `HEAL`.
If you are confident that your mesh is defined correctly (e.g. checked once with the integrity check)
you can disable this check (via `DISABLE`) to avoid the additional runtime overhead of the check.

```python
cube_polyhedron = Polyhedron(
  polyhedral_source=(cube_vertices, cube_faces),
  density=cube_density,
  normal_orientation=NormalOrientation.INWARDS, # OUTWARDS (default) or INWARDS
  integrity_check=PolyhedronIntegrity.VERIFY,   # VERIFY (default), DISABLE or HEAL
)
```

> [!TIP]
> More examples and plots are depicted in the
[jupyter notebook](script/polyhedral-gravity.ipynb).


### Minimal C++ Example

The following example shows how to use the C++ library to compute the gravity.
It works analogously to the Python example above.

```cpp
// Defining the input like above in the Python example
std::vector<std::array<double, 3>> vertices = ...
std::vector<std::array<size_t, 3>> faces = ...
double density = 1.0;
// The constant density polyhedron is defined by its vertices & faces
// It also supports the hand-over of NormalOrientation and PolyhedronIntegrity as optional arguments
// as above described for the Python Interface
Polyhedron polyhedron{vertices, faces, density};
std::vector<std::array<double, 3>> points = ...
std::array<double, 3> point = points[0];
bool parallel = true;
```

The C++ library provides also two ways to compute the gravity. Via
the free function `evaluate`...

```cpp
const auto[pot, acc, tensor] = GravityModel::evaluate(polyhedron, point, parallel);
```

... or via the `GravityEvaluable` class.

```cpp
// Instantiation of the GravityEvaluable object
GravityEvaluable evaluable{polyhedron};

// From now, we can evaluate the gravity model for any point with
const auto[potential, acceleration, tensor] = evaluable(point, parallel);
// or for multiple points with
const auto results = evaluable(points, parallel);
```

Similarly to Python, the C++ implementation also provides mesh checking capabilities.

> [!TIP]
> For reference, have a look at [the main method](./src/main.cpp) of the C++ executable.

## Installation

### With conda

The python interface can be easily installed with
[conda](https://anaconda.org/conda-forge/polyhedral-gravity-model):

```bash
conda install -c conda-forge polyhedral-gravity-model
```

### With pip

As a second option, you can also install the python interface with pip from [PyPi](https://pypi.org/project/polyhedral-gravity/).

```bash
pip install polyhedral-gravity
```

Binaries for the most common platforms are available on PyPI including
Windows, Linux and macOS. For macOS and Linux, binaries for
`x86_64` and `aarch64` are provided.
In case `pip` uses the source distribution, please make sure that
you have a C++17 capable compiler and CMake installed.

### From source

The project uses the following dependencies,
all of them are **automatically** set-up via CMake:

- GoogleTest (1.13.0 or compatible), only required for testing
- spdlog (1.13.0 or compatible), required for logging
- tetgen (1.6 or compatible), required for I/O
- yaml-cpp (0.8.0 or compatible), required for I/O
- thrust (2.1.0 or compatible), required for parallelization and utility
- xsimd (11.1.0 or compatible), required for vectorization of the `atan(..)`
- pybind11 (2.12.0 or compatible), required for the Python interface, but not the C++ standalone

The module will be build using a C++17 capable compiler,
CMake. Just execute the following command in
the repository root folder:

```bash
pip install .
```

To modify the build options (like parallelization) have a look
at the [next paragraph](#building-the-c-library--executable). The options
are modified by setting the environment variables before executing
the `pip install .` command, e.g.:

```bash
export POLYHEDRAL_GRAVITY_PARALLELIZATION="TBB"
pip install .
```

(Optional: For a faster build you can install all dependencies available
for your system in your local python environment. That way, they
won't be fetched from GitHub.)

## C++ Library & Executable

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

The recommended CMake settings using the `TBB` backend would look like this:

```bash
cmake .. -POLYHEDRAL_GRAVITY_PARALLELIZATION="TBB"
```

### Running the C++ Executable

After the build, the gravity model can be run by executing:

```bash
./polyhedralGravity <YAML-Configuration-File>
```

where the YAML-Configuration-File contains the required parameters.
Examples for Configuration Files and Polyhedral Source Files can be
found in this repository in the folder `/example-config/`.

#### Input Configuration File

The configuration should look similar to the given example below.
It is required to specify the source-files of the polyhedron's mesh (more info
about the supported file in the [documentation](https://esa.github.io/polyhedral-gravity-model/quickstart/supported_input.html)), the density
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
    density: 2670.0                             # constant density, units must match with the mesh (see section below)
    points: # Location of the computation point(s) P
      - [ 0, 0, 0 ]                             # Here it is situated at the origin
    check_mesh: true                            # Fully optional, enables mesh autodetect+repair of 
                                                # the polyhedron's vertex ordering (not given: true)
  output:
    filename: "gravity_result.csv"              # The name of the output file 

````

#### Output

The executable produces a CSV file containing $V$, $V_x$, $V_y$, $V_z$,
$V_{xx}$, $V_{yy}$, $V_{zz}$, $V_{xy}$, $V_{xz}$, $V_{yz}$ for every
computation point *P*.

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
