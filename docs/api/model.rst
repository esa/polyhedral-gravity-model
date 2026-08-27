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


Polyhedral Mesh
---------------

The vertices and the triangular faces of a :cpp:class:`polyhedralGravity::Polyhedron` are held by a
:cpp:class:`polyhedralGravity::PolyhedralMesh` as `Kokkos <https://kokkos.org/>`__ views in the memory of
the compute device. Every view is declared with its exact extents, i.e. only the number of vertices or
faces is a runtime dimension while the trailing dimensions are compile time constants, and its layout is
that of a C-contiguous array. A mesh can therefore be built directly on top of the buffer of NumPy,
PyTorch, or JAX, no matter whether that buffer lives on the host or on an accelerator
(see :cpp:enum:`polyhedralGravity::MemoryLocation`).

.. doxygenclass:: polyhedralGravity::PolyhedralMesh

.. doxygenenum:: polyhedralGravity::MemoryLocation


Kokkos Backend
--------------

The namespace :cpp:any:`polyhedralGravity::kokkos` contains the standalone Kokkos parts of this library:
the runtime session (``util/KokkosSession.h``) and the mesh views (``model/PolyhedralMeshView.h``). The
kernels themselves live next to the model code they implement, i.e. the per-face kernel in
``model/GravityModelDetail.h`` and the kernel launches in ``model/GravityEvaluable.cpp``. A user of this
library never has to interact with any of it, since the Kokkos runtime is initialized on demand and
finalized when the process exits.

The mesh views come as a small hierarchy: a :cpp:struct:`polyhedralGravity::kokkos::PolyhedralMeshView`
holds the elementary properties of a polyhedron, i.e. its vertices and faces, and belongs to the
:cpp:class:`polyhedralGravity::Polyhedron`. A
:cpp:struct:`polyhedralGravity::kokkos::GravitationalMeshView` extends it by the caches of Tsoulis'
algorithm, which only make sense inside a :cpp:class:`polyhedralGravity::GravityEvaluable` and are
therefore allocated there — in the memory space of the compute backend that
:cpp:class:`polyhedralGravity::GravityEvaluable` was created for, and only there.

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