.. _examples-python:

Examples Python
===============

Details about mesh and input units can be found in :ref:`quick-start-io`.

The use of the Python interface is pretty straight-forward since
there is only one method: :code:`evaluate(..)` or the alternative
class :code:`GravityEvaluable` caching the polyhedron.
Have a look at :ref:`evaluable-vs-eval` for further
details about the difference.
If you strive for maximal performance, use the the class :code:`GravityEvaluable`.
If you quickly want to try out something, use :code:`evaluate(..)`.
The polyhedral source can either be a tuple of vertices and faces, or
a list of polyhedral mesh files (see :ref:`supported-polyhedron-source-files`).

The method calls follow the same pattern as the C++ interface.

.. code-block:: python

    import polyhedral_gravity as model

    #############################################
    # Free function call
    #############################################
    # Define the polyhedral_source which is either (vertices, faces) or a list of files
    # Define the density (a float)
    # Define the computation_points which is either a single point or a list of points
    # Define if the computation is parallel or not (the default is parallel, which corresponds to True)
    results = model.evaluate(
        polyhedral_source=polyhedral_source,
        density=density,
        computation_points=computation_points,
        parallel=True
    )

    #############################################
    # With the evaluable class
    # This allows multiple evaluations on the same polyhedron
    # without the overhead of setting up the normals etc.
    #############################################
    # Parameters are the same as for the free function call
    evaluable = model.GravityEvaluable(
        polyhedral_source=polyhedral_source,
        density=density
    )
    results = evaluable(
        computation_points=computation_points,
        parallel=True
    )


Free function
-------------

**Example 1:** Evaluating the gravity model for a given polyhedron
defined from within source code for a specific point and density.

.. code-block:: python

        import polyhedral_gravity as model

        # Defining every input parameter in the source code
        vertices = ...          # (N-3)-array-like of type float
        faces = ...             # (N-3)-array-like of type int
        density = ...           # float
        computation_point = ... # (3)-array-like

        # Evaluate the gravity model
        # Notice that the third argument could also be a list of points
        # Returns a tuple of potential, acceleration and tensor
        # If computation_point would be a (N,3)-array, the output would be list of triplets!
        potential, acceleration, tensor = model.evaluate((vertices, faces), density, computation_point, parallel=True)


**Example 2a:** Evaluating the gravity model for a given polyhedron
in some source files for several points and density.
Of course, you can also use keyword arguments for the parameters.

.. code-block:: python

        import polyhedral_gravity as model

        # Reading the vertices and files from a .node and .face file
        file_vertices = '___.node'      # str, path to file
        file_nodes = '___.face'         # str, path to file
        density = ...                   # float
        computation_points = ...        # (N,3)-array-like

        # Evaluate the gravity model
        # Notice that the last argument could also be a list of points
        # Returns a list of tuple of potential, acceleration and tensor
        results = model.evaluate(
            polyhedral_source=[file_vertices, file_nodes],
            density=density,
            computation_points=computation_points,
            parallel=True
        )


**Example 2b:** Evaluating the gravity model for a given polyhedron
in some source files for a specific point and density.

.. code-block:: python

        import polyhedral_gravity as model

        # Reading the vertices and files from a single .mesh file
        file = '___.mesh'       # str, path to file
        density = ...           # float
        computation_point = ... # (3)-array-like

        # Evaluate the gravity model
        # Notice that the last argument could also be a list of points
        # Returns a tuple of potential, acceleration and tensor
        # If computation_point would be a (N,3)-array, the output would be list of triplets!
        potential, acceleration, tensor = model.evaluate([mesh], density, computation_point)


For example 2a and 2b, refer to :ref:`supported-polyhedron-source-files` to view the available
options for polyhedral input.

**Example 3:** A guard statement checks that the plane unit
normals are pointing outwards and no triangular surface is degenerated.
Only use this statement if one needs clarification
about the vertices' ordering due to its quadratic complexity!

.. code-block:: python

    import polyhedral_gravity as model
    import polyhedral_gravity.utility as mesh_sanity

    # Defining every input parameter in the source code
    vertices = ...          # (N-3)-array-like of type float
    faces = ...             # (N-3)-array-like of type int
    density = ...           # float
    computation_point = ... # (3)-array-like



    # Additional guard statement to check that the plane normals
    # are outwards pointing
    if mesh_sanity.check_mesh(vertices, faces):
        # Evaluate the gravity model
        # Returns a tuple of potential, acceleration and tensor
        # If computation_point would be a (N,3)-array, the output would be list of triplets!
        potential, acceleration, tensor = model.evaluate((vertices, faces), density, computation_point)


GravityEvaluable
----------------

Use the :code:`GravityEvaluable` class to cache the polyhedron data over multiple calls.
This drastically improves the performance, as the polyhedral data is "stored" on the C++ side,
rather than being converted from Python to C++ for every call.

This approach is especially useful one wants to calculate multiple points for the same polyhedron, but
the points are not known in advance (e.g. when propagating a spacecraft).
Have a look at the example below to see how to use the :code:`GravityEvaluable` class.

.. code-block:: python

        import polyhedral_gravity as model

        # Defining every input parameter in the source code
        vertices = ...           # (N-3)-array-like of type float
        faces = ...              # (N-3)-array-like of type int
        density = ...            # float
        computation_points = ... # (N,3)-array-like

        # Create the evaluable object
        evaluable = model.GravityEvaluable(polyhedral_source, density)

        for point in computation_points:
            # Evaluate the gravity model for single points (3)-array-like
            potential, acceleration, tensor = evaluable(point, parallel=True)

        # Due to the GravityEvaluable's caching the above for-loop is nearly
        # as fast as the following (find the runtime details below), which returns
        # a list of triplets comprising potential, acceleration, tensor
        results = evaluable(computation_points, parallel=True)
