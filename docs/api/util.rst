Util
====

Overview
--------

The namespace :code:`polyhedralGravity::util` contains utility
for operations on iterable Containers and Constants.

The overloads taking a :code:`std::array` are annotated with
:code:`KOKKOS_INLINE_FUNCTION`, so the gravity model's math runs unchanged
inside a Kokkos kernel, i.e. also on the GPU.

Documentation
-------------

.. doxygennamespace:: polyhedralGravity::util
