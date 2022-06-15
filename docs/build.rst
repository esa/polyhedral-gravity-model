Requirements and Build
======================

CMake
-----

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

============================================== ===================================================================================================================================
Name (Default)                                 Options
============================================== ===================================================================================================================================
PARALLELIZATION_HOST (:code:`CPP`)             :code:`CPP` = Serial Execution on Host/ :code:`OMP`/ :code:`TBB`  = Parallel Execution on Host with OpenMP or Intel's TBB
PARALLELIZATION_DEVICE (:code:`CPP`)           :code:`CPP` = Serial Execution on Device/ :code:`OMP`/ :code:`TBB`  = Parallel Execution on Device with OpenMP or Intel's TBB
LOGGING_LEVEL (:code:`2`)                      :code:`0` = TRACE/ :code:`1` = DEBUG/ :code:`2` = INFO / :code:`3` = WARN/ :code:`4` = ERROR/ :code:`5` = CRITICAL/ :code:`6` = OFF
USE_LOCAL_TBB (:code:`OFF`)                    Use a local installation of :code:`TBB` instead of setting it up via :code:`CMake`
BUILD_POLYHEDRAL_GRAVITY_DOCS (:code:`OFF`)    Build this documentation
BUILD_POLYHEDRAL_GRAVITY_TESTS (:code:`ON`)    Build the Tests
BUILD_POLYHEDRAL_PYTHON_INTERFACE (:code:`ON`) Build the Python interface
============================================== ===================================================================================================================================


Requirements
------------

Requirements (all of them are automatically set-up via :code:`CMake`):

- GoogleTest 1.11.0 (only required for testing)
- spdlog 1.9.2 (required for logging)
- tetgen 1.6 (required for I/O)
- yaml-cpp 0.7.0 (required for I/O)
- thrust 1.16.0 (required for parallelization and utility)
- xsimd 8.1.0 (required for vectorization of the `atan(..)`)



