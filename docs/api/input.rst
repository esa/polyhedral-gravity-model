Input
=====

Overview
--------

The Polyhedral Gravity Model comes with multiple ways of processing
input. In general, there are the two abstract interfaces
:code:`ConfigSource` and :code:`DataSource` and their
respective implementations :code:`YAMLConfigReader` and
:code:`TetgenAdapter`. The first duet is used for
getting the information what execrably to calculate, whereas
the second duet is used for processing the polyhedral source.

Abstract Interfaces
-------------------

.. doxygenclass:: polyhedralGravity::ConfigSource

.. doxygenclass:: polyhedralGravity::DataSource

Implementations
---------------

.. doxygenclass:: polyhedralGravity::YAMLConfigReader

.. doxygenclass:: polyhedralGravity::TetgenAdapter