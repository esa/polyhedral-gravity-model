Requirements, Build, Installation
=================================

Requirements
------------

Requirements (all of them are automatically set-up via :code:`CMake`):

- GoogleTest 1.13.0 (only required for testing)
- spdlog 1.11.0 (required for logging)
- tetgen 1.6 (required for I/O)
- yaml-cpp 0.7.0 (required for I/O)
- thrust 2.1.0 (required for parallelization and utility)
- xsimd 11.1.0 (required for vectorization of the `atan(..)`)
- pybind11 2.10.4 (required for the Python interface, but not the C++ standalone)


Build with CMake
----------------

The C++ implementation relies on :code:`CMake` as build system.
The requirements (see below) are set-up automatically during
the build process. Use the instructions below to build the project, from the
repository's root directory:

.. code-block::

    mkdir build
    cd build
    cmake .. <options>
    cmake --build .

The available options are the following:

================================================ ===================================================================================================================================
Name (Default)                                   Options
================================================ ===================================================================================================================================
POLYHEDRAL_GRAVITY_PARALLELIZATION (:code:`CPP`) :code:`CPP` = Serial Execution / :code:`OMP` or :code:`TBB`  = Parallel Execution with OpenMP or Intel's TBB
LOGGING_LEVEL (:code:`2`)                        :code:`0` = TRACE/ :code:`1` = DEBUG/ :code:`2` = INFO / :code:`3` = WARN/ :code:`4` = ERROR/ :code:`5` = CRITICAL/ :code:`6` = OFF
USE_LOCAL_TBB (:code:`OFF`)                      Use a local installation of :code:`TBB` instead of setting it up via :code:`CMake`
BUILD_POLYHEDRAL_GRAVITY_DOCS (:code:`OFF`)      Build this documentation
BUILD_POLYHEDRAL_GRAVITY_TESTS (:code:`ON`)      Build the Tests
BUILD_POLYHEDRAL_PYTHON_INTERFACE (:code:`ON`)   Build the Python interface
================================================ ===================================================================================================================================

Build & Installation with pip
-----------------------------

Use pip to install the python interface in your local python runtime.
The module will be build using CMake. Just execute in repository root:

.. code-block::

    pip install .

To modify the build options (like parallelization) have a look
at the :code:`setupy.py`. As simple example, to modify the parallelization,
just set the environment variable like below:

.. code-block::

    export POLYHEDRAL_GRAVITY_PARALLELIZATION="TBB"

Installation with conda
-----------------------

The python interface can be easily installed with `conda <https://anaconda.org/conda-forge/polyhedral-gravity-model>`__:

.. code-block::

    conda install -c conda-forge polyhedral-gravity-model

