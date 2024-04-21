.. _mesh-integrity-check:

Mesh Integrity & Explanation
============================

The polyhedral gravity model by *Tsoulis et al.* is defined using
a polyhedron with :code:`OUTWARDS` pointing plane unit normals.
The results are negated if the plane unit normals are :code:`INWARDS` pointing.

This property correlates to the sorting order of the vertices in the faces.
In the context of most programs, one uses here the terminology *clockwise*
and *anti-clockwise* ordering of vertices in the definition of the faces.
The classification also depends on the coordinate system's definition (right- or left-handed).

To be more illustrative and independent of the concrete coordinate system, we stick
to the orientation of the normals: :code:`OUTWARDS` and :code:`INWARDS` pointing.

For example, given you have the following cube:

.. code-block:: python

    import numpy as np

    CUBE_VERTICES = np.array(
      [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
       [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]
    )
    CUBE_FACES_OUTWARDS = np.array(
      [[1, 3, 2], [0, 3, 1], [0, 1, 5], [0, 5, 4], [0, 7, 3], [0, 4, 7],
       [1, 2, 6], [1, 6, 5], [2, 3, 6], [3, 7, 6], [4, 5, 6], [4, 6, 7]]
    )
    cube_inwards = Polyhedron(
      polyhedral_source=(CUBE_VERTICES, CUBE_FACES_OUTWARDS),
      normal_orientation=NormalOrientation.OUTWARDS,
    )

Its plane unit normals are pointing :code:`OUTWARDS`. However, if you simply invert the order like
in the following, they will point :code:`INWARDS`:

.. code-block:: python

    CUBE_FACES_INWARDS = CUBE_FACES_OUTWARDS[:, [1, 0, 2]]
    cube_outwards = Polyhedron(
      polyhedral_source=(CUBE_VERTICES, CUBE_FACES_OUTWARDS),
      normal_orientation=NormalOrientation.INWARDS,
    )

And then, there is still the possibility that one can define a polyhedron by hand
with mixed orientations! This would lead to entirely different (wrong) results.

So, what can we do?

The :code:`Polyhedron` class checks the pointing direction of every single
normal. This way, the `Polyhedron` ensures correct results even if a mistake occurs during the definition.

Since we are in 3D and might have concave and convex polyhedrons,
the viable option is the `Möller–Trumbore intersection algorithm <https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm>`__.
It checks the amount of intersections the plane unit normal has with the polyhedron.
If its an even number, the normals is :code:`OUTWARDS` pointing, otherwise :code:`INWARDS`.
In the *current implementation*, we implement a naiv version
which takes :math:`O(n^2)` operations - which can get quite expensive for polyhedrons with many faces.

To make this as straightforward to use as possible, we provide four options
for the construction of a polyhedron in the form of the enum :code:`PolyhedronIntegrity`:

+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------+----------------+
| :code:`PolyhedronIntegrity` | Meaning                                                                                                                                                                                                                 | Invalid Polyhedron Possible| Overhead       |
+=============================+=========================================================================================================================================================================================================================+============================+================+
| :code:`AUTOMATIC` (default) | Throws an Exception if the :code:`NormalOrientation` is different or inconsistent than specified AND Prints a short-version of the above explanation to :code:`stdout` informing the user about :math:`O(n^2)` runtime  | NO                         | :math:`O(n^2)` |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------+----------------+
| :code:`VERIFY`              | Throws an Exception if the :code:`NormalOrientation` is different or inconsistent than specified                                                                                                                        | NO                         | :math:`O(n^2)` |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------+----------------+
| :code:`HEAL`                | Modifies the specified :code:`NormalOrientation` and the faces array to a consistent polyhedron if they are inconsistent                                                                                                | NO                         | :math:`O(n^2)` |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------+----------------+
| :code:`DISABLE`             | Disables all checks                                                                                                                                                                                                     | YES                        | None           |
+-----------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------+----------------+
