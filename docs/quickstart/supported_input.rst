.. _supported-polyhedron-source-files:

Supported Polyhedron Source Files
=================================

The implementation supports multiple common mesh formats for
the polyhedral source. These include:

====================== ==================================================== ==================================================================================================================================================
File Suffix            Name                                                 Comment
====================== ==================================================== ==================================================================================================================================================
  `.node` and `.face`                     TetGen's files                     These two files need to be given as a pair to the input. `Documentation of TetGen's files <https://wias-berlin.de/software/tetgen/fformats.html>`__
        `.mesh`                         Medit's mesh files                   Single file containing every needed mesh information.
        `.ply`          The Polygon File format/ Stanfoard Triangle format   Single file containing every needed mesh information. Blender File Format.
        `.off`                          Object File Format                   Single file containing every needed mesh information.
        `.stl`                       Stereolithography format                Single file containing every needed mesh information. Blender File Format.
    `.obj` or `.tab`                   Wvaefront OBJ format                  Single file containing every needed mesh information.
====================== ==================================================== ==================================================================================================================================================

**Notice!** Only the ASCII versions of those respective files are supported! This is especially
important for e.g. the `.ply` files which also can be in binary format.

Good tools to convert your Polyhedron to a supported format (also for interchanging
ASCII and binary format) are e.g.:

- `Meshio <https://github.com/nschloe/meshio>`__ for Python
- `OpenMesh <https://openmesh-python.readthedocs.io/en/latest/readwrite.html>`__ for Python

You can select the metric unit of the mesh to be either in
meter :math:`[m]`, kilometer :math:`[km]`, or unitless :math:`[1]` when constructing
an object of type :code:`Polyhedron`. Have a look at the API or at :ref:`examples-python` and :ref:`examples-cpp`.
In case, the mesh is unitless (and specified as such), the implementation won't multiply the results
with the Gravitational Constant :math:`G` and return them raw.