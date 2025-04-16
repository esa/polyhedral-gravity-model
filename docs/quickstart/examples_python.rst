.. _examples-python:

Examples Python
===============

Details about mesh and input units can be found in :ref:`quick-start-io`.

The use of the Python interface is pretty straight-forward since
there is only one method: :code:`evaluate(..)` or the alternative
class :code:`GravityEvaluable` caching the :code:`Polyhedron` and intermediate results.
Have a look at :ref:`evaluable-vs-eval` for further
details about the difference.
If you strive for maximal performance, use the the class :code:`GravityEvaluable`.
If you quickly want to try out something, use :code:`evaluate(..)`.
The polyhedral source can either be a tuple of vertices and faces, or
a list of polyhedral mesh files (see :ref:`supported-polyhedron-source-files`).

The method calls follow the same pattern as the C++ interface.
You always define a :code:`Polyhedron` and then evaluate the polyhedral
gravity model using either the :code:`GravityEvaluable` or :code:`evaluate(..)`.


Free function
-------------

**Example 1:** Evaluating the gravity model for a given constant density
polyhedron defined from within source code for a single point.

.. code-block:: python

        from polyhedral_gravity import Polyhedron, evaluate

        # Defining every input parameter in the source code
        vertices = ...          # (N-3)-array-like of type float
        faces = ...             # (N-3)-array-like of type int
        density = ...           # float
        computation_point = ... # (3)-array-like

        # Create a constant density polyhedron & Evaluate the gravity model
        # Notice that the third argument could also be a list of points
        # Returns a tuple of potential, acceleration and tensor
        # If computation_point would be a (N,3)-array, the output would be list of triplets!
        polyhedron = Polyhedron((vertices, faces), density)
        potential, acceleration, tensor = evaluate(polyhedron, computation_point, parallel=True)


**Example 2a:** Evaluating the gravity model for a given constant density polyhedron
in some source files for several points.
Of course, you can also use keyword arguments for the parameters.
We also check that the polyhedron's plane unit normals are actually
outwards pointing. This will raise a ValueError if at least on
plane unit normal is inwards pointing.

.. code-block:: python

    from polyhedral_gravity import Polyhedron, evaluate, PolyhedronIntegrity, NormalOrientation

    # Reading the vertices and files from a .node and .face file
    file_vertices = '___.node'      # str, path to file
    file_nodes = '___.face'         # str, path to file
    density = ...                   # float
    computation_points = ...        # (N,3)-array-like

    # Evaluate the gravity model
    # Notice that the last argument could also be a list of points
    # Returns a list of tuple of potential, acceleration and tensor
    polyhedron = Polyhedron(
        polyhedral_source=[file_vertices, file_nodes],
        density=density,
        normal_orientation=NormalOrientation.OUTWARDS,
        integrity_check=PolyhedronIntegrity.VERIFY,
    )
    results = evaluate(
        polyhedron=polyhedron,
        computation_points=computation_points,
        parallel=True,
    )


**Example 2b:** Evaluating the gravity model for a given constant density polyhedron
in some source files for a specific point.
We also check that the polyhedron's plane unit normals are actually
outwards pointing. We don't specify this here, as :code:`OUTWARDS` is the default.
It will raise a ValueError if at least on plane unit normal is inwards pointing.

.. code-block:: python

    from polyhedral_gravity import Polyhedron, evaluate, PolyhedronIntegrity

    # Reading the vertices and files from a single .mesh file
    file = '___.mesh'       # str, path to file
    density = ...           # float
    computation_point = ... # (3)-array-like

    # Evaluate the gravity model
    # Notice that the last argument could also be a list of points
    # Returns a tuple of potential, acceleration and tensor
    # If computation_point would be a (N,3)-array, the output would be list of triplets!
    polyhedron = Polyhedron(
        polyhedral_source=[file],
        density=density,
        integrity_check=PolyhedronIntegrity.VERIFY,
    )
    potential, acceleration, tensor = evaluate(polyhedron, computation_point)


For example 2a and 2b, refer to :ref:`supported-polyhedron-source-files` to view the available
options for polyhedral input.

**Example 3a:** Here explicitly disable the security check.
We **won't get an exception** if the plane unit normals are not
oriented as specified, **but we also don't pay for the check with quadratic runtime complexity!**

.. code-block:: python

    from polyhedral_gravity import Polyhedron, evaluate, PolyhedronIntegrity, NormalOrientation

    # Defining every input parameter in the source code
    vertices = ...          # (N-3)-array-like of type float
    faces = ...             # (N-3)-array-like of type int
    density = ...           # float
    computation_point = ... # (3)-array-like

    # Evaluate the gravity model
    # Returns a tuple of potential, acceleration and tensor
    # If computation_point would be a (N,3)-array, the output would be list of triplets!
    polyhedron = Polyhedron(
        polyhedral_source=(vertices, faces),
        density=density,
        normal_orientation=NormalOrientation.OUTWARDS,
        integrity_check=PolyhedronIntegrity.DISABLE,
    )
    potential, acceleration, tensor = evaluate(polyhedron, computation_point)


**Example 3b:** Here we use the :code:`HEAL` option.
This guarantees a valid polyhedron. But the ordering of the faces array and
the normal_orientation might differ.
And we also need to pay the additional quadratic runtime for the checking algorithmus.

.. code-block:: python

    from polyhedral_gravity import Polyhedron, evaluate, PolyhedronIntegrity, NormalOrientation

    # Defining every input parameter in the source code
    vertices = ...          # (N-3)-array-like of type float
    faces = ...             # (N-3)-array-like of type int
    density = ...           # float
    computation_point = ... # (3)-array-like

    # Actually, the normal_orientation doesn't matter! We could the argument
    # as HEAL guarantees a valid polyhedron
    # but the polyhedron might different properties polyhedron.faces
    # and polyhedron.normal_orientation than specified
    polyhedron = Polyhedron(
        polyhedral_source=(vertices, faces),
        density=density,
        normal_orientation=NormalOrientation.OUTWARDS,
        integrity_check=PolyhedronIntegrity.HEAL,
    )
    potential, acceleration, tensor = evaluate(polyhedron, computation_point)


**Example 4:** Here we focus on physical properties of the mesh and density.
The evaluation's output depends on the input units.
If the density is given in :math:`kg/m^3`, and the mesh is in :math:`m`, then the
potential is :math:`m^2/s^2`.
You can change the length unit in case your mesh input is scaled differently
by using the switch :code:`metric_unit`. It allows the values :code:`METER`,
:code:`KILOMETER`, or :code:`UNITLESS`.
In the last case, the gravitational constant :math:`G` is not multiplied to the output.

.. code-block:: python

    from polyhedral_gravity import Polyhedron, evaluate, MetricUnit

    # Defining every input parameter in the source code
    file = '___.obj'        # str, path to file
    density = ...           # float
    computation_point = ... # (3)-array-like

    # Evaluate gravity model with a polyhedron defined whose mesh is in kilometers
    polyhedron = Polyhedron(
        polyhedral_source=[file],        # Mesh in km
        density=density,                 # Density now in kg/km^3
        metric_unit=MetricUnit.KILOMETER # Alternative: METER (default) or UNITLESS
    )
    potential, acceleration, tensor = evaluate(polyhedron, computation_point)


GravityEvaluable
----------------

Use the :code:`GravityEvaluable` class to cache the polyhedron data over multiple calls.
This drastically improves the performance, as the polyhedral data is "stored" on the C++ side,
rather than being converted from Python to C++ for every call.

This approach is especially useful one wants to calculate multiple points for the same polyhedron, but
the points are not known in advance (e.g. when propagating a spacecraft).
Have a look at the example below to see how to use the :code:`GravityEvaluable` class.

.. code-block:: python

        from polyhedral_gravity import Polyhedron, GravityEvaluable, evaluate, PolyhedronIntegrity

        # Defining every input parameter in the source code
        vertices = ...           # (N-3)-array-like of type float
        faces = ...              # (N-3)-array-like of type int
        density = ...            # float
        computation_points = ... # (N,3)-array-like

        # Definition of the Polyhedron in previous examples
        polyhedron = Polyhedron(
            polyhedral_source=(vertices, faces),
            density=density,
            integrity_check=PolyhedronIntegrity.HEAL,
        )

        # Create the evaluable object
        evaluable = GravityEvaluable(polyhedron, density)

        for point in computation_points:
            # Evaluate the gravity model for single points (3)-array-like
            potential, acceleration, tensor = evaluable(point, parallel=True)

        # Due to the GravityEvaluable's caching the above for-loop is nearly
        # as fast as the following (find the runtime details below), which returns
        # a list of triplets comprising potential, acceleration, tensor
        results = evaluable(computation_points, parallel=True)
