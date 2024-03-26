Requirements, Build, Installation
=================================

Requirements
------------

Requirements (all of them are automatically set-up via :code:`CMake`):

- GoogleTest (1.13.0 or compatible), only required for testing
- spdlog (1.11.0 or compatible), required for logging
- tetgen (1.6 or compatible), required for I/O
- yaml-cpp (0.7.0 or compatible), required for I/O
- thrust (2.1.0 or compatible), required for parallelization and utility
- xsimd (11.1.0 or compatible), required for vectorization of the :code:`atan(..)`
- pybind11 (2.10.4 or compatible), required for the Python interface, but not the C++ standalone


Build with CMake
----------------

The C++ implementation relies on :code:`CMake` as build system.
The requirements (see below) are set-up automatically during
the build process. Use the instructions below to build the project, from the
repository's root directory:

.. code-block:: bash

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

Installation with conda
-----------------------

The python interface can be easily installed with `conda <https://anaconda.org/conda-forge/polyhedral-gravity-model>`__:

.. code-block:: bash

    conda install -c conda-forge polyhedral-gravity-model

Installation with pip
---------------------

As a second option, you can also install the python interface with pip from
`PyPi <https://pypi.org/project/polyhedral-gravity/>`__:

.. code-block:: bash

    pip install polyhedral-gravity

Binaries for the most common platforms are available on PyPI including
Windows, Linux and macOS. For macOS and Linux, binaries for
:code:`x86_64` and :code:`aarch64` are provided.
In case :code:`pip` uses the source distribution, please make sure that
you have a C++17 capable compiler and CMake installed.


Build & Installation from source
--------------------------------

Use pip to install the python interface in your local python runtime.
The module will be build using CMake. Just execute in repository root:

.. code-block::

    pip install .

To modify the build options (like parallelization) have a look
at the :code:`setupy.py`. The options are modified by setting the
environment variables before executing the :code:`pip install .` command, e.g.:

.. code-block::

    export POLYHEDRAL_GRAVITY_PARALLELIZATION="TBB"



