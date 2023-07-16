Quick Start C++
===============

As Executable
-------------

General
~~~~~~~

After the build, the gravity model can be run by executing:

.. code-block::

    ./polyhedralGravity <YAML-Configuration-File>

where the YAML-Configuration-File contains the required parameters.
Examples for Configuration Files and Polyhedral Source Files can be
found in this repository in the folder `/example-config/`.

Config File
~~~~~~~~~~~

The configuration should look similar to the given example below.
It is required to specify the source-files of the polyhedron's mesh (more info
about the supported files below), the density
of the polyhedron, and the wished computation points where the
gravity tensor shall be computed.
Further one must specify the name of the .csv output file.

.. code-block:: yaml

    ---
    gravityModel:
      input:
        polyhedron: #polyhedron source-file(s)
          - "../example-config/data/tsoulis.node"   # .node contains the vertices
          - "../example-config/data/tsoulis.face"   # .face contains the triangular faces
        density: 2670.0                             # constant density in [kg/m^3]
        points: # Location of the computation point(s) P
          - [ 0, 0, 0 ]                             # Here it is situated at the origin
        check_mesh: true                            # Fully optional, enables input checking (not given: false)
      output:
        filename: "gravity_result.csv"              # The name of the output file


Have a look at :ref:`supported-polyhedron-source-files` to view the available
options for polyhedral input.

As Library
----------

The following paragraph gives some examples on how to
use the polyhedral model as library from within source code
(All examples with :code:`using namespace polyhedralGravity`).

Free function
~~~~~~~~~~~~~


**Example 1:** Evaluating the gravity model for a given polyhedron
defined from within source code for a specific point and density.
Further, we disable the parallelization using the optional fourth parameter (which defaults to true).

.. code-block:: cpp

        // Defining every input parameter in the source code
        std::vector<std::array<double, 3>> vertices = ...
        std::vector<std::array<size_t, 3>> faces = ...
        Polyhedron polyhedron{vertices, faces};
        double density = ...
        std::array<double, 3> point = ...

        // Returns either a single of vector of results
        // Here, we only have one point. Thus we get a single result
        const auto[pot, acc, tensor] = GravityModel::evaluate(polyhedron, density, point, false);


**Example 2:** Evaluating the gravity model for a given polyhedron
in some source files for a specific point and density.
Further, we explicitly enable the parallelization using the optional fourth parameter
(which defaults to true).

.. code-block:: cpp

        // Reading just the polyhedral source from file,
        // whereas the rest is defined within source-code
        auto files = std::vector<std::string>{"tsoulis.node", "tsoulis.face"};
        double density = ...
        std::array<double, 3> point = ...

        // Returns either a single of vector of results
        // Here, we only have one point. Thus we get a single result
        const auto[pot, acc, tensor] = GravityModel::evaluate(files, density, point, true);


**Example 3:** Evaluating the gravity model for a given configuration
from a .yaml file.

.. code-block:: cpp

        // Reading the configuration from a yaml file
        std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>("config.yaml");
        Polyhedron poly = config->getDataSource()->getPolyhedron();
        double density = config->getDensity();
        // This time, we use multiple points
        std::vector<std::array<double, 3>> points = config->getPointsOfInterest();

        // Returns either a single of vector of results
        // Here, we have multiple point. Thus we get a vector of results!
        const results = GravityModel::evaluate(poly, density, points);

**Example 4:** A guard statement checks that the plane unit
normals are pointing outwards and no triangle is degenerated.
Only use this statement if one needs clarification
about the vertices' ordering due to its quadratic complexity!

.. code-block:: cpp

        // Reading the configuration from a yaml file
        std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>("config.yaml");
        Polyhedron poly = config->getDataSource()->getPolyhedron();
        double density = config->getDensity();
        std::array<double, 3> point = config->getPointsOfInterest()[0];

        // Guard statement
        if (MeshChecking::checkTrianglesNotDegenerated(poly) && MeshChecking::checkNormalsOutwardPointing(poly)) {
            GravityResult result = GravityModel::evaluate(poly, density, point);
        }



Class-Based method
~~~~~~~~~~~~~~~~~~

Evaluating the gravity model for a given polyhedron
defined from within source code for a specific point and density.

.. code-block:: cpp

        // Defining every input parameter in the source code
        std::vector<std::array<double, 3>> vertices = ...
        std::vector<std::array<size_t, 3>> faces = ...
        Polyhedron polyhedron{vertices, faces};
        double density = ...
        std::array<double, 3> point = ...
        std::vector<std::array<double, 3>> points = ...

        // Instantiation of the GravityEvaluable object
        GravityEvaluable evaluable{polyhedron, density};

        // From now, we can evaluate the gravity model for any point with
        const auto[pot, acc, tensor] = evaluable(point);
        // or for multiple points with
        const auto results = evaluable(points);
        // and we can also disable e.g. the parallelization like for the free function
        const auto singleResultTuple = evaluable(point, false);
