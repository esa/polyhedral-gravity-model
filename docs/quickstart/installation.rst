Installation
============



Installation with conda
-----------------------

The python interface can be easily installed with `conda <https://anaconda.org/conda-forge/polyhedral-gravity-model>`__:

.. code-block:: bash

    conda install -c conda-forge polyhedral-gravity-model

The python package on conda is parallelized with OpenMP.
It is currently available for all operating systems (macOS, Linux, Windows), but
only for :code:`x86_64` systems.
Have a look at the :ref:`installation-pip`. It also provides wheels for :code:`aarch64`.


.. _installation-pip:

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


Installation & Build from source
--------------------------------

All these steps require a working C/C++ Compiler and CMake to be installed.

Build Python Package
~~~~~~~~~~~~~~~~~~~~

Use pip to install the python interface in your local python runtime.
The module will be build using CMake. Just execute in repository root:

.. code-block::

    pip install .

To modify the build options (like parallelization) have a look
at the :ref:`build-options` for an overview of options.
The pip installation will call CMake. To modify build options, just set them as
environment variable before executing the :code:`pip install .` command, e.g.:

.. code-block::

    export POLYHEDRAL_GRAVITY_PARALLELIZATION="TBB"


Build C++ library/ executable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The C++ implementation relies on :code:`CMake` as build system.
The requirements (see below) are set-up automatically during
the build process. Use the instructions below to build the project, from the
repository's root directory:

.. code-block:: bash

    mkdir build
    cd build
    cmake .. <options>
    cmake --build .

Have a look at :ref:`build-options` for an overview of options available for the CMake build.


.. _build-options:

Build Options
~~~~~~~~~~~~~

The available options are the following:

====================================================== ============================================================================================================
Name (Default)                                         Options
====================================================== ============================================================================================================
POLYHEDRAL_GRAVITY_PARALLELIZATION (:code:`CPP`)       :code:`CPP` = Serial Execution / :code:`OMP` or :code:`TBB`  = Parallel Execution with OpenMP or Intel's TBB
POLYHEDRAL_GRAVITY_LOGGING_LEVEL (:code:`INFO`)        :code:`TRACE`, :code:`DEBUG`, :code:`INFO`, :code:`WARN`, :code:`ERROR`, :code:`CRITICAL`, :code:`OFF`
BUILD_POLYHEDRAL_GRAVITY_DOCS (:code:`OFF`)            Build this documentation
BUILD_POLYHEDRAL_GRAVITY_TESTS (:code:`ON`)            Build the Tests
BUILD_POLYHEDRAL_GRAVITY_PYTHON_INTERFACE (:code:`ON`) Build the Python interface
====================================================== ============================================================================================================

Dependencies (automatically set-up)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dependencies (all of them are automatically set-up via :code:`CMake`):

- GoogleTest (1.15.2 or compatible), only required for testing
- spdlog (1.13.0 or compatible), required for logging
- tetgen (1.6 or compatible), required for I/O
- yaml-cpp (0.8.0 or compatible), required for I/O
- thrust (2.1.0 or compatible), required for parallelization and utility
- xsimd (11.1.0 or compatible), required for vectorization of the :code:`atan(..)`
- pybind11 (2.12.0 or compatible), required for the Python interface, but not the C++ standalone

Build this documentation
------------------------

In order to build this documentation from source, you require the following dependencies:

- :code:`Doxygen`
- :code:`Sphinx` with the following plugins
    - :code:`breathe`
    - :code:`sphinx-book-theme`
- The :code:`polyhedral_gravity` Python Package needs to be installed

How you install the :code:`polyhedral_gravity` Python Package is stated above.
The other dependencies can be install them with your favorite package manager (e.g. conda, pip, brew, apt,...):

.. code-block:: bash

    conda install doxygen sphinx breathe sphinx-book-theme


Build the documentation via CMake
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can build the documentation locally using CMake by executing the following commands:

.. code-block:: bash

    mkdir build && cd build
    cmake .. -DBUILD_POLYHEDRAL_GRAVITY_DOCS=ON
    cmake --build . --target Doxygen
    cmake --build . --target Sphinx
    open docs/sphinx/index.html

If you installed the Sphinx dependencies in a non-standard-path, e.g., a conda environment, you might need to help CMake
finding it by specifying during the configure step `-DCMAKE_PREFIX_PATH=$CONDA_PREFIX`

Build the documentation invoking Sphinx
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also omit CMake and build the documentation locally in the following way:

.. code-block:: bash

    cd docs
    export BUILD_DOCS_CLI=1
    make html # Alternatively: sphinx-build -M html . ./_build
    open _build/html/index.html
