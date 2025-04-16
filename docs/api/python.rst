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

    Specifies the parallelization of the :code:`polyhedral_gravity` module at compile time
    (One of :code:`CPP`, :code:`OMP`, :code:`TBB`)

    This attribute corresponds to the value of the ``POLYHEDRAL_GRAVITY_PARALLELIZATION`` C++ variable.

.. py:attribute:: __commit__

    Specifies the Git commit hash of the source code at compile time.

    It is assigned the value of the ``POLYHEDRAL_GRAVITY_COMMIT_HASH`` C++ variable.

.. py:attribute:: __logging__

    Specifies the logging level fixed at compile time.

    This corresponds to the value defined by the ``POLYHEDRAL_GRAVITY_LOGGING_LEVEL`` C++ variable.
