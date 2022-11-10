Quick Start Python
==================

The use of the Python interface is pretty straight-forward since
there is only one method: :code:`evaluate(..)`.

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
        # Notice that the last argument could also be a list of points
        # Returns (list of) tuple of potential, acceleration and tensor
        potential, acceleration, tensor = model.evaluate(vertices, faces, density, computation_point)


**Example 2a:** Evaluating the gravity model for a given polyhedron
in some source files for a specific point and density.

.. code-block:: python

        import polyhedral_gravity as model

        # Reading the vertices and files from a .node and .face file
        file_vertices = '___.node'      # str, path to file
        file_nodes = '___.face'         # str, path to file
        density = ...                   # float
        computation_point = ...         # [] or np.array of length 3

        # Evaluate the gravity model
        # Notice that the last argument could also be a list of points
        # Returns (list of) tuple of potential, acceleration and tensor
        potential, acceleration, tensor = model.evaluate([file_vertices, file_nodes], density, computation_point)


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
        # Returns (list of) tuple of potential, acceleration and tensor
        potential, acceleration, tensor = model.evaluate([mesh], density, computation_point)


For example 2a and 2b, refer to :ref:`supported-polyhedron-source-files` to view the available
options for polyhedral input.

**Example 3:** A guard statement checks that the plane unit
normals are pointing outwards. Only use this statement if one needs clarification
about the vertices' ordering due to its quadratic complexity!

.. code-block:: python

    import polyhedral_gravity as model
    import polyhedral_gravity.utility as model_sanity

    # Defining every input parameter in the source code
    vertices = ...              # [] of [] or np.array of length 3 and type float
    faces = ...                 # [] of [] or np.array of length 3 and type int
    density = ...               # float
    computation_point = ...     # [] or np.array of length 3

    # Evaluate the gravity model
    # Notice that the last argument could also be a list of points
    # Returns (list of) tuple of potential, acceleration and tensor

    # Additional guard statement
    if model_sanity.check_input(vertices, faces):
        potential, acceleration, tensor = model.evaluate(vertices, faces, density, computation_point)

