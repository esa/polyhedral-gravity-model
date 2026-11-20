Input
=====

Overview
--------

The Polyhedral Gravity Model comes with multiple ways of processing
input. In general, there are is one abstract interface
:code:`ConfigSource` to specify a model evaluation.
Its implementation is realized in the class :code:`YAMLConfigReader`.

The file input is delegated from the :code:`YAMLConfigReader` to
the functions in the namespace :code:`MeshReader` which
vice-versa delegates calls to Tetgen's file formats to the
:code:`TetgenAdapter`.


Configuration Input
-------------------

.. doxygenclass:: polyhedralGravity::ConfigSource

.. doxygenclass:: polyhedralGravity::YAMLConfigReader

Mesh File Input
---------------

.. doxygennamespace:: polyhedralGravity::MeshReader

.. doxygenclass:: polyhedralGravity::TetgenAdapter