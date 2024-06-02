Implementation's Details
========================


This code is a validated implementation in C++17 of the Polyhedral Gravity Model
by Tsoulis et al.. Additionally, the model provides a Python binding.
It was initially created in a collaborative project between
TU Munich and ESA's Advanced Concepts Team.

If this implementation proves useful to you, please consider citing the
`accompanying paper <https://doi.org/10.21105/joss.06384>`__
published in the *Journal of Open Source Software*.
It summarizes application scenarios, related work and the architectural design since version :code:`3.0`.

For a detailed dive into the theoretical background, consult the
`initial project report <https://mediatum.ub.tum.de/doc/1695208/1695208.pdf>`__.
Its content about the theoretical background, core implementation, verification, and results is still valid.
However, the implementation's public interfaces changed over time.
The initial project report refers to version :code:`1.2.1`.
Hence, its usage examples will not work with the current version!

Further, the approach is documented by Tsoulis et al. in the following two major publications:

- `Tsoulis, D., 2012. Analytical computation of the full gravity tensor of a homogeneous arbitrarily shaped polyhedral source using line integrals. Geophysics, 77(2), pp.F1-F11. <http://dx.doi.org/10.1190/geo2010-0334.1>`__ and its corresponding implementation in FORTRAN_ (last accessed: 12.09.2022).
- `TSOULIS, Dimitrios; GAVRIILIDOU, Georgia. A computational review of the line integral analytical formulation of the polyhedral gravity signal. Geophysical Prospecting, 2021, 69. Jg., Nr. 8-9, S. 1745-1760 <https://doi.org/10.1111/1365-2478.13134>`__ and its corresponding implementation in MATLAB_ (last accessed: 28.03.2024), which is strongly based on the former implementation in FORTRAN.



.. _FORTRAN: https://software.seg.org/2012/0001/index.html
.. _MATLAB: https://github.com/Gavriilidou/GPolyhedron