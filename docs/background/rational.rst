Statement of Need
=================

The Rational
------------

.. image:: ../../_static/eros_010.png
  :align: center
  :width: 400
  :alt: Downscaled mesh of (433) Eros to 10% of its original vertices and faces.

Downscaled mesh of (433) Eros to 10% of its original vertices and faces.

The complex gravitational fields of irregular bodies, such as asteroids and comets,
are often modeled using polyhedral gravity models since alternative approaches like
mascon models or spherical harmonics struggle with these bodies' irregular geometry.
The spherical harmonics approach struggles with convergence close to the surface [1]_,
whereas mascon models require a computationally expensive amount of mascons
(point masses of which the target body comprises) to model fine-granular surface geometry [2]_.

In contrast, polyhedral gravity models provide an analytic solution for the computation of the
gravitational potential, acceleration (and second derivative) given
a mesh of the body [3]_ [4]_ with the only assumption of homogeneous density.
The computation of the gravitational potential and acceleration is a computationally expensive task,
especially for large meshes, which can however benefit from parallelization either over computed target
points for which we seek potential and acceleration or over the mesh. Thus, a high-performance
implementation of a polyhedral gravity model is desirable.

While some research code for these models exists, they are not focused on usability and are
limited to FORTRAN_ and proprietary software like MATLAB_.
There is a lack of well-documented, actively maintained open-source implementations,
particularly in modern programming languages, and with a focus on scalability and performance.

This circumstance and the fact that polyhedral models are often used in studying gravitational
fields, e.g., for Eros [5]_, or as a reference for creating
new neural models [6]_ make an easy-to-install implementation necessary.

The presented software has already seen application in several research works.
It has been used to optimize trajectories around the highly irregular comet 67P/Churyumov-Gerasimenko
with the goal of maximizing the gravity signal [7]_ using pygmo [8]_.
In the context of that work, the presented implementation was extended to enable caching and even
serialization to persistent memory on the C++ side. A change that enables researchers to, e.g.,
efficiently propagate an orbit since the computation points can be given apiece and do not need
to be all known from the beginning.

Further, it has been used to study the effectiveness of so-called
neural density fields [9]_, where it served as ground truth to
(pre-)train neural networks representing the density
distribution of an arbitrarily shaped body [10]_.

Thus, this model is highly versatile overall due to its easy-to-use API.
It can be used in a wide range of applications, especially due to the availability
on major platforms like Windows, macOS, and Linux for ARM64 and x86_64.
We hope it will enable further research in the field, especially related to
recent machine-learning techniques, which typically rely on Python implementations.


References
----------

.. [1] Šprlák, M., & Han, S.-C. (2021). On the use of spherical harmonic series inside the minimum brillouin sphere: Theoretical review and evaluation by GRAIL and LOLA satellite data. Earth-Science Reviews, 222, 103739. https://doi.org/10.1016/j.earscirev.2021.103739
.. [2] Wittick, P. T., & Russell, R. P. (2017). Mascon models for small body gravity fields. AAS/AIAA Astrodynamics Specialist Conference, 162, 17–162.
.. [3] Tsoulis, D. (2012). Analytical computation of the full gravity tensor of a homogeneous arbitrarily shaped polyhedral source using line integrals. Geophysics, 77(2), F1–F11. https://doi.org/10.1190/geo2010- 0334.1
.. [4] Tsoulis, D., & Gavriilidou, G. (2021). A computational review of the line integral analytical formulation of the polyhedral gravity signal. Geophysical Prospecting, 69(8-9), 1745–1760. https://doi.org/10.1111/1365- 2478.13134
.. [5] Zhang, Z., Cui, H., Cui, P., & Yu, M. (2010). Modeling and analysis of gravity field of 433Eros using polyhedron model method. 2010 2nd International Conference on Information Engineering and Computer Science, 1–4. https://doi.org/10.1109/iciecs.2010.5677738
.. [6] Martin, J., & Schaub, H. (2023). The physics-informed neural network gravity model revisited: Model generation III. 33rd AAS/AIAA Space Flight Mechanics Meeting, Austin, United States.
.. [7] Maråk, R., Blazquez, E., & Gómez, P. (2023). Trajectory optimization of a spacecraft swarm orbiting around 67P/Churyumov-Gerasimenko. Proceedings of the 9th International Conference on Astrodynamics Tools and Techniques, ICATT. https://doi.org/10.5270/ esa-gnc-icatt-2023-057
.. [8] Biscani, F., & Izzo, D. (2020). A parallel global multiobjective framework for optimization: pagmo. Journal of Open Source Software, 5(53), 2338. https://doi.org/10.21105/joss.02338
.. [9] Izzo, D., & Gómez, P. (2022). Geodesy of irregular small bodies via neural density fields. Communications Engineering, 1(1), 48. https://doi.org/10.1038/s44172-022-00050-3
.. [10] Schuhmacher, J., Gratl, F., Izzo, D., & Gómez, P. (2023). Investigation of the robustness of neural density fields. Proceedings of the 12th International Conference on Guidance, Navigation & Control Systems (GNC). https://doi.org/10.5270/esa-gnc-icatt-2023-067


.. _FORTRAN: https://software.seg.org/2012/0001/index.html
.. _MATLAB: https://github.com/Gavriilidou/GPolyhedron