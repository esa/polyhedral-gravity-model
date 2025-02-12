.. _examples-cpp:

Examples C++
============

Details about mesh and input units can be found in :ref:`quick-start-io`.

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
        density: 2670.0                             # constant density, units must match with the mesh (see section below)
                                                    # Depends on metric_unit: 'km' -> kg/km^3, 'm' -> kg/m^3, 'unitless' -> 'unitless'
        points: # Location of the computation point(s) P
          - [ 0, 0, 0 ]                             # Here it is situated at the origin
        check_mesh: true                            # Fully optional, enables mesh autodetect+repair of
                                                    # the polyhedron's vertex ordering (not given: true)
        metric_unit: m                              # One of 'm', 'km' or ' unitless' (not given: 'm')
      output:
        filename: "gravity_result.csv"              # The name of the output file


Have a look at :ref:`supported-polyhedron-source-files` to view the available
options for polyhedral input.


As Library
----------

The following paragraph gives some examples on how to
use the polyhedral model as library from within source code
(All examples with :code:`using namespace polyhedralGravity`).

Individual Function (without caching)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


**Example 1:** Evaluating the gravity model for a given polyhedron
defined from within source code for a specific point and density.
Further, we disable the parallelization using the optional fourth parameter (which defaults to true)
and we check the if the plane unit normals are actually outwards pointing
(costing :math:`O(n^2)`).

.. code-block:: cpp

        // Defining every input parameter in the source code
        std::vector<std::array<double, 3>> vertices = ...
        std::vector<std::array<size_t, 3>> faces = ...
        double density = ...
        std::array<double, 3> point = ...

        // Returns either a single of vector of results
        // Here, we only have one point. Thus we get a single result
        Polyhedron polyhedron{vertices, faces, density, NormalOrientation::OUTWARDS, PolyhedronIntegrity::VERIFY};
        const auto[pot, acc, tensor] = GravityModel::evaluate(polyhedron, point, false);


**Example 2:** Evaluating the gravity model for a given polyhedron
in some source files for a specific point and density.
Further, we explicitly enable the parallelization using the optional fourth parameter
(which defaults to true). We disable the :math:`O(n^2)` check if the
plane unit normals are actually outwards pointing.

.. code-block:: cpp

        // Reading just the polyhedral source from file,
        // whereas the rest is defined within source-code
        auto files = std::vector<std::string>{"tsoulis.node", "tsoulis.face"};
        double density = ...
        std::array<double, 3> point = ...

        // Returns either a single of vector of results
        // Here, we only have one point. Thus we get a single result
        Polyhedron polyhedron{files, density, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE};
        const auto[pot, acc, tensor] = GravityModel::evaluate(polyhedron, point, true);


**Example 3:** Evaluating the gravity model for a given configuration
from a .yaml file.

.. code-block:: cpp

        // Reading the configuration from a yaml file
        std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>("config.yaml");
        auto polyhedralSource = config->getDataSource()->getPolyhedralSource();
        double density = config->getDensity();
        PolyhedronIntegrity checkPolyhedralInput = config->getMeshInputCheckStatus() ? PolyhedronIntegrity::HEAL : PolyhedronIntegrity::DISABLE;
        // This time, we use multiple points
        std::vector<std::array<double, 3>> points = config->getPointsOfInterest();

        // Returns either a single of vector of results
        // Here, we have multiple point. Thus we get a vector of results!
        Polyhedron polyhedron{polyhedralSource, density, NormalOrientation::OUTWARDS, checkPolyhedralInput};
        const results = GravityModel::evaluate(polyhedron, points);

**Example 4:** If our :code:`Polyhedron` contains any inconsistencies in the definition, e.g.,
some plane unit normals point outwards and some point inwards,
the :code:`HEAL` option will fix the :code:`NormalOrientation` and any inconsistencies
related to it. It won't throw an exception.
The result will always be fine.
Further, you can also specify the metric unit of the mesh, which is one of :code:`MetricUnit`,
either :code:`METER`, :code:`KILOMETER`, or :code:`UNITLESS`.
Important! This also influences the unit of the output and the required density input.

.. code-block:: cpp

        // Reading the configuration from a yaml file
        std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>("config.yaml");
        auto polyhedralSource = config->getDataSource()->getPolyhedralSource();
        double density = config->getDensity();
        std::array<double, 3> point = config->getPointsOfInterest()[0];
        const auto metricUnit = config->getMeshUnit();

        Polyhedron polyhedron{polyhedralSource, density, NormalOrientation::OUTWARDS, PolyhedronIntegrity::HEAL, metricUnit};
        const auto[pot, acc, tensor] = GravityModel::evaluate(polyhedron, point);


GravityEvaluable (with caching)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Evaluating the gravity model for a given polyhedron
defined from within source code for a specific point and density.


.. code-block:: cpp

        // Defining every input parameter in the source code
        std::vector<std::array<double, 3>> vertices = ...
        std::vector<std::array<size_t, 3>> faces = ...
        Polyhedron polyhedron{vertices, faces};
        double density = ...
        Polyhedron polyhedron{vertices, faces, density, NormalOrientation::OUTWARDS, PolyhedronIntegrity::DISABLE};

        // Our computation points
        std::array<double, 3> point = ...
        std::vector<std::array<double, 3>> points = ...

        // Instantiation of the GravityEvaluable object
        GravityEvaluable evaluable{polyhedron};

        // From now, we can evaluate the gravity model for any point with
        const auto[pot, acc, tensor] = evaluable(point);
        // or for multiple points with
        const auto results = evaluable(points);
        // and we can also disable e.g. the parallelization like for the free function
        const auto singleResultTuple = evaluable(point, false);
