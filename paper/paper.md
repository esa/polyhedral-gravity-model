---
title: 'Efficient Polyhedral Gravity Modeling in Modern C++ and Python'
tags:
  - C++
  - Python
  - astronomy
  - dynamics
  - asteroids
  - gravity
  - numerical methods
  - polyhedral model
authors:
  - name: Jonas Schuhmacher
    orcid: 0009-0005-9693-4530
    affiliation: 1
  - name: Emmanuel Blazquez
    orcid: 0000-0001-9697-582X
    affiliation: 2
  - name: Fabio Gratl
    orcid: 0000-0001-5195-7919
    affiliation: 1
  - name: Dario Izzo
    orcid: 0000-0002-9846-8423
    affiliation: 2
  - name: Pablo Gómez
    orcid: 0000-0002-5631-8240
    affiliation: 2
affiliations:
 - name: Chair for Scientific Computing, Technische Universität München, Arcisstraße 21, 80333 München, Germany 
   index: 1
 - name: Advanced Concepts Team, European Space Agency, European Space Research and Technology Centre (ESTEC), Keplerlaan 1, 2201 AZ Noordwijk, The Netherlands
   index: 2
date: 07 July 2023
bibliography: paper.bib

---

# Summary

Polyhedral gravity models are ubiquitous for modelling the gravitational field of irregular bodies, such as asteroids and comets. We present an open-source C++ library for the efficient, parallelized computation of a polyhedral gravity model. We also provide a Python interface to the library using *pybind11*. The library is particularly focused on delivering high performance and scalability, which we achieve through vectorization and parallelization with *xsimd* and *thrust*, respectively. The library supports many common formats, such as *.stl*, *.off*, *.ply*, *.mesh* and *tetgen*'s *.node* and *.face*.

# Statement of Need

The complex gravitational fields of irregular bodies, such as asteroids and comets, are often modeled using polyhedral gravity models as they provide an analytic solution for the computation of the gravitational potential, acceleration (and second derivative) given a mesh of the body [@tsoulis2012analytical;@tsoulis2021computational]. The computation of the gravitational potential and acceleration is a computationally expensive task, especially for large meshes, which can however benefit from parallelization either over computed target points for which we seek potential and acceleration or over the mesh. Thus, a high-performance implementation of a polyhedral gravity model is desirable. 

While some research code for these models exists, they are not focused on usability and are limited to FORTRAN \footnote{\url{https://software.seg.org/2012/0001/index.html}} and proprietary software like MATLAB \footnote{\url{https://github.com/Gavriilidou/GPolyhedron}}. There is a lack of well-documented, actively maintained open-source implementations, particularly in modern programming languages, and with a focus on scalability and performance.

The presented software has already seen application in several research works. It has been used to optimize trajectories around the highly irregular comet 67P/Churyumov-Gerasimenko [@marak2023trajectory]. Further, it has been used to study the effectiveness of so-called neural density fields [@izzo2022geodesy], where it can serve as ground truth and to pre-train neural networks [@GeodesyNetsBenchmark]. **TODO_add_more_examples**

Thus, overall this model is highly versatile and can be used in a wide range of applications. We hope it will enable further research in the field, especially related to recent machine learning techniques, which typically rely on Python implementations.

# Polyhedral Model

On a mathematical level, the implemented model follows the approach by Petrović [@petrovic1996determination] as refined by Tsoulis and Petrović [@tsoulis2001singularities]. A comprehensive description of the mathematical foundations of the model is given in the associated student report [@jonas_report]. 

Implementation-wise, it makes use of the inherent parallelization opportunity of the approach as it iterates over the mesh. This parallelization is achieved via *thrust*, which allows utilizing *OpenMP* and *Intel TBB*. On a finer scale, individual, costly operations were investigated and, e.g., the \texttt{arctan} operations were vectorized using *xsimd*. On an application side, the user may use the implemented caching mechanism to avoid the recomputation of mesh properties, such as the face normals.

Extensive tests using GoogleTest are used via GitHub Actions to ensure the (continued) correctness of the implementation.

# Installation \& Contribution

The library is available on GitHub \footnote{\url{https://github.com/esa/polyhedral-gravity-model}} and can be installed with *pip* or from *conda* \footnote{\url{https://anaconda.org/conda-forge/polyhedral-gravity-model}}. Build instructions using *CMake* are provided in the repository. The library is licensed under a GPL license.

The project is open to contributions via pull requests with instructions on how to contribute provided in the repository.

# Usage Instructions

We provide detailed usage instructions in the technical documentation on ReadTheDocs \footnote{\url{https://polyhedral-gravity-model-cpp.readthedocs.io/en/latest/}}. Additionally, a minimal working example is given in the repository readme and more extensive examples as a *Jupyter* notebook \footnote{\url{https://github.com/esa/polyhedral-gravity-model/blob/main/script/polyhedral-gravity.ipynb}}.

# References
