# polyhedral-gravity-model

![Build and Test](https://github.com/schuhmaj/polyhedral-gravity-model-cpp/actions/workflows/ctest.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/polyhedral-gravity-model-cpp/badge/?version=latest)](https://polyhedral-gravity-model-cpp.readthedocs.io/en/latest/?badge=latest)

Implementation of the Polyhedral Gravity Model in C++ 17.

The implementation is based on the paper [Tsoulis, D., 2012. Analytical computation of the full gravity tensor of a homogeneous arbitrarily shaped polyhedral source using line integrals. Geophysics, 77(2), pp.F1-F11.](http://dx.doi.org/10.1190/geo2010-0334.1)
and its corresponding [implementation in FORTRAN](https://software.seg.org/2012/0001/index.html).

Supplementary details can be found in the more recent paper [TSOULIS, Dimitrios; GAVRIILIDOU, Georgia. A computational review of the line integral analytical formulation of the polyhedral gravity signal. Geophysical Prospecting, 2021, 69. Jg., Nr. 8-9, S. 1745-1760.](https://doi.org/10.1111/1365-2478.13134)
and its corresponding [implementation in MATLAB](https://github.com/Gavriilidou/GPolyhedron),
which is strongly based on the former implementation in FORTRAN.

## Documentation

The full extensive documentation can be found on [readthedocs](https://polyhedral-gravity-model-cpp.readthedocs.io/en/latest/).


## Requirements
The project uses the following dependencies, 
all of them are **automatically** set-up via CMake:

- GoogleTest 1.11.0 (only required for testing)
- spdlog 1.9.2 (required for logging)
- tetgen 1.6 (required for I/O)
- yaml-cpp 0.7.0 (required for I/O)
- thrust 1.16.0 (required for parallelization and utility)
- xsimd 8.1.0 (required for vectorization of the `atan(..)`)

## Python interface

Use pip to install the python interface in your local python runtime.
The module will be build using CMake and the using the above
requirements. Just execute in repository root:

    pip install .

To modify the build options (like parallelization) have a look
at the `setupy.py` and the [next paragraph](#build-c).

## Build C++

The program is build by using CMake. So first make sure that you installed
CMake and then follow these steps:

    mkdir build
    cd build
    cmake .. <options>
    cmake --build .

The following options are available:

|         Name (Default)         |                                                       Options and Remark                                                      |
|:------------------------------:|:-----------------------------------------------------------------------------------------------------------------------------:|
|  PARALLELIZATION_HOST (`CPP`)  |             `CPP` = Serial Execution on Host, `OMP`/ `TBB` = Parallel Execution on Host with OpenMP or Intel's TBB              |
| PARALLELIZATION_DEVICE (`CPP`) | `CPP`= Serial Execution on Device, `OMP`/ `TBB` = Parallel Execution on Device with OpenMP or Intel's TBB |
|      LOGGING_LEVEL (`2`)       |                          `0`= TRACE, `1`=DEBUG, `2`=INFO, `3`=WARN, `4`=ERROR, `5`=CRITICAL, `6`=OFF                          |

During testing the combination PARALLELIZATION_HOST=`CPP` and PARALLELIZATION_DEVICE=`TBB`
has been the most performant.
It is further not recommend to change the LOGGING_LEVEL to something else than `INFO=2`.

The recommended CMake command would look like this (we only need to change `PARALLELIZATION_DEVICE`, since
the defaults of the others are already correctly set):

    cmake .. -DPARALLELIZATION_DEVICE="TBB"

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
    polyhedron:                                 #polyhedron source-file(s)
      - "../example-config/data/tsoulis.node"   #.node contains the vertices
      - "../example-config/data/tsoulis.face"   #.face contains the triangular faces
    density: 2670.0                             #constant density in [kg/m^3]
    points:                                     #Location of the computation point(s) P
      - [0, 0, 0]                               #Here it is situated at the origin
    check: true                                 #Fully optional, input checking is enabled
  output:
    filename: "gravity_result.csv"              #The name of the output file 

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

The polyhedron's plane unit normals must point outwards!
Typically, the nodes of a mesh are either given in clockwise or counter-clockwise ordering.
To check that you've the right format the program is equipped with an input checking procedure using the
Möller–Trumbore algorithm. However, this procedure has quadratic complexity. Hence, it can be optionally enabled
the yaml configuration file.

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
