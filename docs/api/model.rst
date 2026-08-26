Model
=====

Overview
--------

The model is the heart of the Polyhedral Gravity Model
since it contains the the two major classes:

* :cpp:class:`polyhedralGravity::Polyhedron`
* :cpp:class:`polyhedralGravity::GravityEvaluable`

The class :cpp:class:`polyhedralGravity::Polyhedron` stores the mesh and density
information and governs the compliance with Tsoulis et al.’s gravity model’s preconditions
It ensures that all plane unit normals of a constructed Polyhedron are consistently
pointing :cpp:enumerator:`polyhedralGravity::NormalOrientation::OUTWARDS` or
:cpp:enumerator:`polyhedralGravity::NormalOrientation::INWARDS`.
Depending on the value of :cpp:enum:`polyhedralGravity::PolyhedronIntegrity`, these
constraints are either enforced by modifying input mesh data or by throwing
an :cpp:class:`std::invalid_argument` exception.
For this purpose, it uses the `Möller–Trumbore intersection algorithm <https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm>`__.


The class :cpp:class:`polyhedralGravity::GravityEvaluable` is provides as a way to
perform the evaluation of the polyhedral gravity model repeatedly
without the need to re-initialize the polyhedron and the gravity model as
caching is performed.
It takes a :cpp:class:`polyhedralGravity::Polyhedron`, uploads it to the compute devices, and provides an
:cpp:func:`polyhedralGravity::GravityEvaluable::operator()` to evaluate the
model at computation point(s) :math:`P` on the
:cpp:enum:`polyhedralGravity::ComputeBackend` of one's choice.
The evaluation itself is implemented once with `Kokkos <https://kokkos.org/>`__ and therefore runs
unchanged on a single CPU thread, on all CPU threads, and on a GPU. Its floating point precision is
selected with :cpp:enum:`polyhedralGravity::ComputePrecision` when constructing the evaluable.
A :cpp:func:`polyhedralGravity::GravityModel::evaluate` summarizes the
functionality :cpp:class:`polyhedralGravity::GravityEvaluable`, but does not
provide any caching throughout multiple calls.


Polyhedron
----------

.. doxygenclass:: polyhedralGravity::Polyhedron

.. doxygenenum:: polyhedralGravity::NormalOrientation

.. doxygenenum:: polyhedralGravity::PolyhedronIntegrity

.. doxygenenum:: polyhedralGravity::MetricUnit


GravityModel
------------

.. doxygenclass:: polyhedralGravity::GravityEvaluable

.. doxygenenum:: polyhedralGravity::ComputeBackend

.. doxygenenum:: polyhedralGravity::ComputePrecision

.. doxygennamespace:: polyhedralGravity::GravityModel


Kokkos Backend
--------------

The namespace :cpp:any:`polyhedralGravity::kokkos` contains everything which directly touches Kokkos:
the runtime session, the device-resident polyhedron, and the kernels. A user of this library never has
to interact with it, since the Kokkos runtime is initialized on demand and finalized when the process
exits.

.. doxygennamespace:: polyhedralGravity::kokkos


Named Tuple
-----------

.. doxygenstruct:: polyhedralGravity::DistanceTemplate

.. doxygenstruct:: polyhedralGravity::TranscendentalExpressionTemplate

.. doxygenstruct:: polyhedralGravity::HessianPlaneTemplate

.. doxygentypedef:: polyhedralGravity::Distance

.. doxygentypedef:: polyhedralGravity::TranscendentalExpression

.. doxygentypedef:: polyhedralGravity::HessianPlane

Type Definitions
----------------

.. doxygentypedef:: polyhedralGravity::Vector3

.. doxygentypedef:: polyhedralGravity::Vector6

.. doxygentypedef:: polyhedralGravity::Vector3Triplet

.. doxygentypedef:: polyhedralGravity::Array3

.. doxygentypedef:: polyhedralGravity::Array6

.. doxygentypedef:: polyhedralGravity::IndexArray3

.. doxygentypedef:: polyhedralGravity::Array3Triplet

.. doxygentypedef:: polyhedralGravity::GravityModelResult

.. doxygentypedef:: polyhedralGravity::PolyhedralFiles

.. doxygentypedef:: polyhedralGravity::PolyhedralSource