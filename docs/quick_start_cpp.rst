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
        check: true                                 # Fully optional, enables input checking (not given: false)
      output:
        filename: "gravity_result.csv"              # The name of the output file


Have a look at :ref:`supported-polyhedron-source-files` to view the available
options for polyhedral input.

As Library
----------

The following paragraph gives some examples on how to
use the polyhedral model as library from within source code
(All examples with :code:`using namespace polyhedralGravity`).


**Example 1:** Evaluating the gravity model for a given polyhedron
defined from within source code for a specific point and density.

.. code-block:: cpp

        // Defining every input parameter in the source code
        std::vector<std::array<double, 3>> vertices = ...
        std::vector<std::array<size_t, 3>> faces = ...
        double density = ...
        std::array<double, 3> point = ...

        // Main method, notice that the last argument can also be a
        // vector of points
        GravityResult result = GravityModel::evaluate({vertices, faces}, density, point);


**Example 2:** Evaluating the gravity model for a given polyhedron
in some source files for a specific point and density.

.. code-block:: cpp

        // Reading just the polyhedral source from file,
        // whereas the rest is defined within source-code
        TetgenAdapter input{{"tsoulis.node", "tsoulis.face"}};
        Polyhedron poly = input.getPolyhedron();
        double density = ...
        std::array<double, 3> point = ...

        // Main method, notice that the last argument can also be a
        // vector of points
        GravityResult result = GravityModel::evaluate(poly, density, point);


**Example 3:** Evaluating the gravity model for a given configuration
from a .yaml file.

.. code-block:: cpp

        // Reading the configuration from a yaml file
        std::shared_ptr<ConfigSource> config = std::make_shared<YAMLConfigReader>("config.yaml");
        Polyhedron poly = config->getDataSource()->getPolyhedron();
        double density = config->getDensity();
        std::array<double, 3> point = config->getPointsOfInterest()[0];

        // Main method, notice that the last argument can also be a
        // vector of points
        GravityResult result = GravityModel::evaluate(poly, density, point);
