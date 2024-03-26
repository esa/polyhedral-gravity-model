.. _evaluable-vs-eval:

GravityEvaluable vs. Evaluate Function
======================================

Below is a comparison of the performance of the standalone free function and the evaluable class.
The benchmark was conducted with a M1 Pro 10-Core CPU (ARM64) using the TBB backend.
The calculation consisted each of 1000 points for the mesh of Eros (24235 vertices and 14744 faces).
The results are shown in the table below (the lower the better/ faster).
Basically, as soon as you have more than one point to evaluate, the evaluable class is faster and
thus recommended. This result comes from the fact that the polyhedral data is cached on the C++ side
and thus does not need to be converted from Python to C++ for every call. Further, the evaluable class
also caches the normals and the volume of the polyhedron, which is not the case for the standalone function.

+----------------------------------------+-------------------------------+
| Test                                   | Time Per Point (microseconds) |
+========================================+===============================+
| Free Function (1000 times 1 point)     | 7765.073                      |
+----------------------------------------+-------------------------------+
| Free Function (1 time 1000 points)     | 275.917                       |
+----------------------------------------+-------------------------------+
| GravityEvaluable (1000 times 1 Point)  | 313.408                       |
+----------------------------------------+-------------------------------+
| GravityEvaluable (1 time 1000 Points)  | 253.031                       |
+----------------------------------------+-------------------------------+

As takeaway, one should always use the :code:`GravityEvaluable` due to its caching
properties. The free function is maily of advantage for few points/ quickly trying
out something.