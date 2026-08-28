Python API - polyhedral_gravity
===============================

Module Overview
---------------

.. automodule:: polyhedral_gravity

Polyhedron
----------

.. autoclass:: polyhedral_gravity.Polyhedron
   :members:
   :special-members: __init__, __getitem__, __repr__

Enums to specify Mesh(-checks)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: polyhedral_gravity.NormalOrientation

.. autoclass:: polyhedral_gravity.PolyhedronIntegrity

.. autoclass:: polyhedral_gravity.MetricUnit


GravityModel
------------

Enums to specify where and how to evaluate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: polyhedral_gravity.ComputeBackend

.. autoclass:: polyhedral_gravity.ComputePrecision


Single Function
~~~~~~~~~~~~~~~

.. autofunction:: polyhedral_gravity.evaluate


Cached Evaluation
~~~~~~~~~~~~~~~~~

.. autoclass:: polyhedral_gravity.GravityEvaluable
   :members:
   :special-members: __init__, __call__, __repr__


Embedded Information
--------------------

The module also provides several attributes that define key metadata
about the build configuration, versioning, and runtime behavior of the
:code:`polyhedral_gravity` module which are set at compile time.
Below is the list of the available attributes:

.. py:attribute:: __version__

    Specifies the version of the :code:`polyhedral_gravity` module at compile time.

    This corresponds to the value of the ``POLYHEDRAL_GRAVITY_VERSION`` C++ variable.

.. py:attribute:: __parallelization__

    Lists the `Kokkos <https://kokkos.org/>`__ execution spaces the :code:`polyhedral_gravity` module was
    compiled with, as a comma-separated string, e.g. :code:`'Serial, OpenMP, Cuda'`.

    These are the backends a :py:class:`polyhedral_gravity.ComputeBackend` can be served by:
    :code:`Serial` for :code:`CPU_SERIAL`, :code:`OpenMP` for :code:`CPU_PARALLEL`, and
    :code:`Cuda`/ :code:`HIP`/ :code:`SYCL` for :code:`GPU_PARALLEL`. If no GPU space is listed,
    requesting :code:`GPU_PARALLEL` raises a :code:`RuntimeError`.

.. py:attribute:: __commit__

    Specifies the Git commit hash of the source code at compile time.

    It is assigned the value of the ``POLYHEDRAL_GRAVITY_COMMIT_HASH`` C++ variable.

.. py:attribute:: __logging__

    Specifies the logging level fixed at compile time.

    This corresponds to the value defined by the ``POLYHEDRAL_GRAVITY_LOGGING_LEVEL`` C++ variable.

.. py:attribute:: __host_compiler__

    Specifies the compiler which translated the host code, e.g. :code:`'GNU 13.2.0'`.

    This corresponds to the value of the ``POLYHEDRAL_GRAVITY_HOST_COMPILER`` C++ variable.

.. py:attribute:: __device_compiler__

    Specifies the compiler which translated the device code, e.g. :code:`'NVIDIA nvcc 12.4'`.
    It is :code:`'None'` if the module was compiled without a GPU backend, i.e. if requesting
    :code:`ComputeBackend.GPU_PARALLEL` raises a :code:`RuntimeError`.

    This corresponds to the value of the ``POLYHEDRAL_GRAVITY_DEVICE_COMPILER`` C++ variable.


PyTorch Interface (Differentiable)
-----------------------------------

The optional :code:`polyhedral_gravity.torch` submodule provides a pure-PyTorch,
autograd-differentiable re-implementation of :code:`evaluate(..)`. It is not
imported by default and requires PyTorch to be installed separately
(:code:`pip install torch`). See :ref:`examples-python` for a usage example,
and the `PyTorch interface notebook <https://github.com/esa/polyhedral-gravity-model/blob/main/examples/notebooks/polyhedral-gravity-torch.ipynb>`__
and `benchmark notebook <https://github.com/esa/polyhedral-gravity-model/blob/main/examples/notebooks/polyhedral-gravity-torch-benchmark.ipynb>`__
for further details.

.. py:function:: polyhedral_gravity.torch.evaluate(vertices, faces, density, computation_points, gravitational_constant=6.67430e-11)

    Gravitational potential, acceleration, and gradient tensor for a polyhedron,
    differentiable w.r.t. vertex positions and density.

    :param vertices: ``(N, 3)`` vertex positions [m] (``torch.Tensor``).
    :param faces: ``(F, 3)`` triangle vertex indices, integer dtype (``torch.Tensor``).
    :param density: Constant density [kg/m^3].
    :param computation_points: ``(Q, 3)`` evaluation positions [m] (``torch.Tensor``).
    :param gravitational_constant: Gravitational constant, defaults to ``6.67430e-11``.
    :returns: Tuple of ``potential (Q,)``, ``acceleration (Q, 3)``, and gradient ``tensor (Q, 6)``
        (ordered as ``[Vxx, Vyy, Vzz, Vxy, Vxz, Vyz]``).
