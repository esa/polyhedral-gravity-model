Calculation
===========

Overview
--------

The calculation is the heart of the Polyhedral Gravity Model
since it contains the important function :code:`evaluate(..)`
used evaluation of the polyhedral gravity model at computation
point(s) P.

Additionally, the class :code:`GravityEvaluable` is provided as a way to
perform the evaluation of the polyhedral gravity model repeatedly
without the need to re-initialize the polyhedron and the gravity model as
some caching is performed. It is recommended to use this feature when evaluating the model repeatedly on a small number of points as the overhead of setting up the models is larger than the computational cost.

Further, it implements the `Möller–Trumbore intersection algorithm <https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm>`__
capable of checking if the plane unit normals of the polyhedron are pointing outwards.

GravityEvaluable
----------------

.. doxygenclass:: polyhedralGravity::GravityEvaluable


GravityModel
------------

.. doxygennamespace:: polyhedralGravity::GravityModel


.. _mesh-checking-cpp:

MeshChecking
------------

.. doxygennamespace:: polyhedralGravity::MeshChecking

Other
-----

.. doxygenfunction:: polyhedralGravity::transformPolyhedron
