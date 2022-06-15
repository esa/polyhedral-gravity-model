Python API
==========

.. py::module:: polyhedral_gravity

Computes the full gravity tensor for a given constant density polyhedron which consists of some vertices and
triangular faces at a given computation point :math:`P`

.. py:function:: evaluate(vertices, faces, density, computation_point)
   :noindex:

   Evaluate the full gravity tensor for a given constant density polyhedron which consists of some vertices and
   triangular faces at a given computation point P

   :param List[List[float[3]]] vertices: vertices of the polyhedron
   :param List[List[int[3]]] faces: faces of the polyhedron
   :param float density: constant density in :math:`\frac{kg}{m^3}`
   :param List[float[3]] computation_point: cartesian computation point :math:`P`
   :return: tuple of potential, acceleration, and second derivative tensor
   :rtype: Tuple[float, List[float[3]], List[float[6]]]


.. py:function:: evaluate(vertices, faces, density, computation_points)
   :noindex:

   Evaluate the full gravity tensor for a given constant density polyhedron which consists of some vertices and
   triangular faces at multiple given computation points

   :param List[List[float[3]]] vertices: vertices of the polyhedron
   :param List[List[int[3]]] faces: faces of the polyhedron
   :param float density: constant density in :math:`\frac{kg}{m^3}`
   :param List[List[float[3]]] computation_points: multiple cartesian computation points :math:`P`
   :return: list of tuple of potential, acceleration, and second derivative tensor
   :rtype: List[Tuple[float, List[float[3]], List[float[6]]]]

.. py:function:: evaluate(input_files, density, computation_point)
   :noindex:

   Evaluate the full gravity tensor for a given constant density polyhedron which consists of some vertices and
   triangular faces at a given computation point :math:`P`

   :param List[str] input_files: polyhedral source files
   :param float density: constant density in :math:`\frac{kg}{m^3}`
   :param List[float[3]] computation_point: cartesian computation point :math:`P`
   :return: tuple of potential, acceleration, and second derivative tensor
   :rtype: Tuple[float, List[float[3]], List[float[6]]]

.. py:function:: evaluate(input_files, density, computation_points)
   :noindex:

    Evaluate the full gravity tensor for a given constant density polyhedron which consists of some vertices and
    triangular faces at multiple given computation points

   :param List[str] input_files: polyhedral source files
   :param List[List[int[3]]] faces: faces of the polyhedron
   :param float density: constant density in :math:`\frac{kg}{m^3}`
   :param List[List[float[3]]] computation_points: multiple cartesian computation points :math:`P`
   :return: list of tuple of potential, acceleration, and second derivative tensor
   :rtype: List[Tuple[float, List[float[3]], List[float[6]]]]
