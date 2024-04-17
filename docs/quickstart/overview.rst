.. _quick-start-io:

Quick Start (C++ & Python)
==========================


This section describes the specification of required inputs and output like
required units and characteristics of the mesh.
These constraints apply for both, the C++ and Python interface.
The actual utilization is presented in :ref:`examples-python` and :ref:`examples-cpp`.


Gravity Model Input
-------------------


The evaluation of the polyhedral gravity model requires the following parameters.

+------------------------------------------------------------------------------+
| Name                                                                         |
+==============================================================================+
| Polyhedral Mesh (either as vertices & faces or as polyhedral source files)   |
+------------------------------------------------------------------------------+
| Constant Density :math:`\rho`                                                |
+------------------------------------------------------------------------------+

The polyhedron's mesh's units must match with the constant density!
For example, if the mesh is in :math:`[m]`, then the constant density should be in :math:`[\frac{kg}{m^3}]`.

Have a look at :ref:`supported-polyhedron-source-files` for an overview of supported polyhedral
mesh files.

.. note::

    The plane unit normals of every face of the polyhedral mesh must point
    consistently **outwards** or **inwards** the polyhedron!
    This property is automatically enforced by the class :code:`Polyhedron` in
    both C++ library and Python interface as long as not explicit set to :code:`DISABLE`.
    Setting this to :code:`DISABLE`, is recommended for advanced users or when you "know your mesh".
    For details, refer to the APIs or :ref:`mesh-integrity-check`.


Gravity Model Output
--------------------

The calculation outputs the following parameters for every Computation Point *P*.
The units of the respective output depend on the units of the input parameters (mesh and density)!

Hence, if e.g. your mesh is in :math:`km`, the density must match. Further, the output units will match the input units.

+------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+-----------------------------------------------------------------+
|         Name                                                                                   | If mesh :math:`[m]` and density :math:`[\frac{kg}{m^3}]`                   |                             Comment                             |
+================================================================================================+============================================================================+=================================================================+
|         :math:`V`                                                                              |  :math:`\frac{m^2}{s^2}` or :math:`\frac{J}{kg}`                           |           The potential or also called specific energy          |
+------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+-----------------------------------------------------------------+
|     :math:`V_x`, :math:`V_y`, :math:`V_z`                                                      |   :math:`\frac{m}{s^2}`                                                    |The gravitational acceleration in the three cartesian directions |
+------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+-----------------------------------------------------------------+
| :math:`V_{xx}`, :math:`V_{yy}`, :math:`V_{zz}`, :math:`V_{xy}`, :math:`V_{xz}`, :math:`V_{yz}` |   :math:`\frac{1}{s^2}`                                                    |The spatial rate of change of the gravitational acceleration     |
+------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+-----------------------------------------------------------------+

This model's output obeys to the geodesy and geophysics sign conventions.
Hence, the potential :math:`V` for a polyhedron with a mass :math:`m > 0` is defined as **positive**.

The accelerations :math:`V_x`, :math:`V_y`, :math:`V_z` are defined as

.. math::

    \textbf{g} = + \nabla V = \left( \frac{\partial V}{\partial x}, \frac{\partial V}{\partial y}, \frac{\partial V}{\partial z} \right)

Accordingly, the second derivative tensor is defined as the derivative of :math:`\textbf{g}`.
