Output
======

Overview
--------

The Polyhedral Gravity Model in general prints the output for
all given computation points to a .csv file in order to keep the terminal clean.
In order to achieve this the :code:`CSVWriter` is used.

This module also contains the :code:`PolyhedralGravityLogger` which serves
as a Wrapper class for accessing :code:`spdlog`'s Logger.

Implementations
---------------

.. doxygenclass:: polyhedralGravity::CSVWriter

.. doxygenclass:: polyhedralGravity::PolyhedralGravityLogger
