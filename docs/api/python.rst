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

Enums to specify MeshChecks
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: polyhedral_gravity.NormalOrientation

.. autoclass:: polyhedral_gravity.PolyhedronIntegrity


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