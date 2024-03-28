Implementation's Details
========================


This code is a validated implementation in C++17 of the Polyhedral Gravity Model
by Tsoulis et al.. It was created in a collaborative project between
TU Munich and ESA's Advanced Concepts Team. Please refer to the
`project report <https://mediatum.ub.tum.de/doc/1695208/1695208.pdf>`__
for extensive information about the theoretical background, related work,
implementation & design decisions, application, verification,
and runtime measurements of the presented code.

The implementation is based on the
paper `Tsoulis, D., 2012. Analytical computation of the full gravity tensor of a homogeneous arbitrarily shaped polyhedral source using line integrals. Geophysics, 77(2), pp.F1-F11. <http://dx.doi.org/10.1190/geo2010-0334.1>`__
and its corresponding implementation in FORTRAN_.

Supplementary details can be found in the more recent
paper `TSOULIS, Dimitrios; GAVRIILIDOU, Georgia. A computational review of the line integral analytical formulation of the polyhedral gravity signal. Geophysical Prospecting, 2021, 69. Jg., Nr. 8-9, S. 1745-1760 <https://doi.org/10.1111/1365-2478.13134>`__
and its corresponding implementation in MATLAB_,
which is strongly based on the former implementation in FORTRAN.



.. _FORTRAN: https://software.seg.org/2012/0001/index.html
.. _MATLAB: https://github.com/Gavriilidou/GPolyhedron