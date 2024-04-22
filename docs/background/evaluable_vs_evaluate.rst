.. _evaluable-vs-eval:

GravityEvaluable vs. Evaluate Function
======================================

The :code:`GravityEvaluable` was introduced with version :code:`v2.0`.
Hence, the `project report <https://mediatum.ub.tum.de/doc/1695208/1695208.pdf>`__ does not contain any further details.
This sections closes the knowledge gap and presents a runtime comparison.

Below is a comparison of the performance of the standalone free function and the evaluable class, as of :code:`v3.0`.
The benchmark was conducted with a M1 Pro 10-Core CPU (ARM64) using the *TBB backend*.
The calculation consisted each of 1000 points for the mesh of Eros (24235 vertices and 14744 faces).
The points are either handed over one-by-one as :math:`(3)`-arrays or as one :math:`(1000, 3)`-array.
The former user case appears when the latter locations are unknown beforehand - like in a trajectory propagation scenario.
The results are shown in the plot below (the lower the better/ faster).

Basically, as soon as you have more than one point to evaluate, the evaluable class is faster and
thus recommended. This result comes from the fact that the polyhedral data is cached on the C++ side
and thus does not need to be converted from Python to C++ for every call. Further, the evaluable class
also caches the normals and the volume of the polyhedron, which is not the case for the standalone function.
As of :code:`v3.0`, the evaluate free function is nearly as good as the :code:`GravityEvaluable` if all points are
known beforehand, since the :code:`Polyhedron` class and its construction are now separated from the
evaluation of the gravity model.

.. image:: /../_static/runtime_measurements.png
  :align: center
  :width: 100%
  :alt: Runtime Measurements for the mesh of (433) Eros in :code:`v3.0`

If the scenario is yet to be determined, use the :code:`GravityEvaluable` due to its caching
properties. If you know all points "beforehand", the approach does not matter!