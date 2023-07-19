Quick Start Python
==================

The use of the Python interface is pretty straight-forward since
there is only one method: :code:`evaluate(..)` or the alternative
class :code:`GravityEvaluable` caching the polyhedron data.

The method calls follow the same pattern as the C++ interface. Thus it is always:

.. code-block:: python

    import polyhedral_gravity as model

    #############################################
    # Free function call
    #############################################
    # Define the polyhedral_source which is either (vertices, faces) or a list of files
    # Define the density (a float)
    # Define the computation_points which is either a single point or a list of points
    # Define if the computation is parallel or not (the default is parallel, which corresponds to True)
    results = model.evaluate(polyhedral_source, density, computation_point, parallel)

    #############################################
    # With the evaluable class
    # This allows multiple evaluations on the same polyhedron
    # without the overhead of setting up the normals etc.
    #############################################
    # Parameters are the same as for the free function call
    evaluable = model.GravityEvaluable(polyhedral_source, density)
    results = evaluable(computation_point, parallel)


Free function
-------------

**Example 1:** Evaluating the gravity model for a given polyhedron
defined from within source code for a specific point and density.

.. code-block:: python

        import polyhedral_gravity as model

        # Defining every input parameter in the source code
        vertices = ...          # [] of [] or np.array of length 3 and type float
        faces = ...             # [] of [] or np.array of length 3 and type int
        density = ...           # float
        computation_point = ... # [] or np.array of length 3

        # Evaluate the gravity model
        # Notice that the third argument could also be a list of points
        # Returns a tuple of potential, acceleration and tensor
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
        computation_points = ...        # [] or np.array of length 3

        # Evaluate the gravity model
        # Notice that the last argument could also be a list of points
        # Returns a list of tuple of potential, acceleration and tensor
        results = model.evaluate(
            polyhedral_source=[file_vertices, file_nodes],
            density=density,
            computation_points=computation_points
            parallel=True
        )


**Example 2b:** Evaluating the gravity model for a given polyhedron
in some source files for a specific point and density.

.. code-block:: python

        import polyhedral_gravity as model

        # Reading the vertices and files from a single .mesh file
        file = '___.mesh'       # str, path to file
        density = ...           # float
        computation_point = ... # [] or np.array of length 3

        # Evaluate the gravity model
        # Notice that the last argument could also be a list of points
        # Returns a tuple of potential, acceleration and tensor
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
    vertices = ...              # [] of [] or np.array of length 3 and type float
    faces = ...                 # [] of [] or np.array of length 3 and type int
    density = ...               # float
    computation_point = ...     # [] or np.array of length 3

    # Evaluate the gravity model
    # Notice that the last argument could also be a list of points
    # Returns (list of) tuple of potential, acceleration and tensor

    # Additional guard statement
    if mesh_sanity.check_mesh(vertices, faces):
        potential, acceleration, tensor = model.evaluate((vertices, faces), density, computation_point)


Evaluable class
---------------

Use the :code:`GravityEvaluable` class to cache the polyhedron data over multiple calls.
This drastically improves the performance, as the polyhedral data is "stored" on the C++ side,
rather than being converted from Python to C++ for every call.

This approach is especially useful one wants to calculate multiple points for the same polyhedron, but
the points are not known in advance (e.g. when propagating a spacecraft).
Have a look at the example below to see how to use the :code:`GravityEvaluable` class.

.. code-block:: python

        import polyhedral_gravity as model

        # Defining every input parameter in the source code
        vertices = ...           # [] of [] or np.array of length 3 and type float
        faces = ...              # [] of [] or np.array of length 3 and type int
        density = ...            # float
        computation_points = ... # [] or np.array of length 3

        # Create the evaluable object
        evaluable = model.GravityEvaluable(polyhedral_source, density)

        for point in computation_points:
            # Evaluate the gravity model
            potential, acceleration, tensor = evaluable(point, parallel=True)

Below is a comparison of the performance of the free function and the evaluable class.
Basically, as soon as you have more than one point to evaluate, the evaluable class is faster and
thus recommended.

+----------------------------------------+-------------------------------+
| Test                                   | Time Per Point (microseconds) |
+========================================+===============================+
| Free Function (1000 times 1 point)     | 7765.073                      |
+----------------------------------------+-------------------------------+
| Free Function (1 time N points)        | 275.917                       |
+----------------------------------------+-------------------------------+
| GravityEvaluable (1000 times 1 Point)  | 313.408                       |
+----------------------------------------+-------------------------------+
| GravityEvaluable (1 time 1000 Points)  | 253.031                       |
+----------------------------------------+-------------------------------+
