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
you have a C++20 capable compiler and CMake installed.


Installation & Build from source
--------------------------------

All these steps require a working C/C++ Compiler and CMake to be installed.

Build Python Package
~~~~~~~~~~~~~~~~~~~~

Use pip to install the python interface in your local python runtime.
The module will be build using CMake. Just execute in repository root:

.. code-block::

    pip install .

The build is driven by `scikit-build-core <https://github.com/scikit-build/scikit-build-core>`__,
which calls CMake. Have a look at :ref:`build-options` for an overview of the options.
Any of them can be set per install with pip's :code:`--config-settings`, abbreviated :code:`-C`:

.. code-block:: bash

    pip install . -C cmake.define.POLYHEDRAL_GRAVITY_DEVICE_BACKEND=CUDA

.. code-block:: bash

    pip install . -C cmake.define.POLYHEDRAL_GRAVITY_FAST_MATH=ON

For :code:`POLYHEDRAL_GRAVITY_DEVICE_BACKEND`, :code:`POLYHEDRAL_GRAVITY_LOGGING_LEVEL`, and
:code:`POLYHEDRAL_GRAVITY_FAST_MATH` the equally named environment variables of the previous,
setuptools-based build keep working as well:

.. code-block:: bash

    export POLYHEDRAL_GRAVITY_DEVICE_BACKEND="CUDA"
    pip install .


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

====================================================== ==========================================================================================================================
Name (Default)                                         Options
====================================================== ==========================================================================================================================
POLYHEDRAL_GRAVITY_DEVICE_BACKEND (:code:`AUTO`)       :code:`AUTO` = detect the GPU vendor's paradigm / :code:`NONE` = CPU-only build / :code:`CUDA`, :code:`HIP`, or :code:`SYCL`
POLYHEDRAL_GRAVITY_LOGGING_LEVEL (:code:`INFO`)        :code:`TRACE`, :code:`DEBUG`, :code:`INFO`, :code:`WARN`, :code:`ERROR`, :code:`CRITICAL`, :code:`OFF`
POLYHEDRAL_GRAVITY_FAST_MATH (:code:`OFF`)             Faster, less accurate FLOAT32 arithmetic. FLOAT64 is bit-identical either way
BUILD_POLYHEDRAL_GRAVITY_DOCS (:code:`OFF`)            Build this documentation
BUILD_POLYHEDRAL_GRAVITY_TESTS (:code:`ON`)            Build the Tests
BUILD_POLYHEDRAL_GRAVITY_PYTHON_INTERFACE (:code:`ON`) Build the Python interface
====================================================== ==========================================================================================================================

The host backends need no configuration: the Kokkos :code:`Serial` backend is always compiled in, and the
:code:`OpenMP` backend is enabled whenever an OpenMP installation is found. On macOS this is Homebrew's
:code:`libomp`; without it, :code:`ComputeBackend.CPU_PARALLEL` falls back to serial execution.

:code:`POLYHEDRAL_GRAVITY_DEVICE_BACKEND=AUTO` picks the GPU vendor's native paradigm from the installed
toolchain: CUDA if :code:`nvcc` is found, HIP if the ROCm compiler is found, SYCL for the Intel LLVM
compiler, and no GPU backend otherwise.

Kokkos compiles the device code with the **C++ compiler**, not with a separate CUDA/ HIP compiler, so
:code:`CXX` has to be one that understands the paradigm: a CUDA-capable :code:`clang++` or NVIDIA's
:code:`nvcc_wrapper` for CUDA, :code:`hipcc` for HIP, and :code:`icpx` for SYCL. If :code:`clang++`
rejects your CUDA version, either point it at an older toolkit with
:code:`-DKokkos_CUDA_DIR=/path/to/cuda` or build with :code:`nvcc_wrapper`:

.. code-block:: bash

    CXX=/path/to/nvcc_wrapper cmake .. -DPOLYHEDRAL_GRAVITY_DEVICE_BACKEND=CUDA

Kokkos determines the GPU architecture by running a probe at configure time, which needs a visible GPU.
On a login node without one, name the architecture yourself instead:

.. code-block:: bash

    cmake .. -DKokkos_ARCH_AMPERE80=ON

Dependencies (automatically set-up)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dependencies (all of them are automatically set-up via :code:`CMake`):

- GoogleTest (1.15.2 or compatible), only required for testing
- spdlog (1.13.0 or compatible), required for logging
- tetgen (1.6 or compatible), required for I/O
- yaml-cpp (0.8.0 or compatible), required for I/O
- Kokkos (5.1.1 or compatible), required for the parallelization on the CPU and the GPU
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
