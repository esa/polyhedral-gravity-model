.. Polyhedral Gravity Model documentation master file, created by
   sphinx-quickstart on Wed May 25 15:10:43 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Polyhedral Gravity Model's documentation!
====================================================

.. image:: https://joss.theoj.org/papers/10.21105/joss.06384/status.svg
   :target: https://doi.org/10.21105/joss.06384
.. image:: https://img.shields.io/github/license/esa/polyhedral-gravity-model
   :alt: GitHub License
.. image:: https://img.shields.io/github/actions/workflow/status/esa/polyhedral-gravity-model/.github%2Fworkflows%2Fbuild-and-test.yml?logo=GitHub%20Actions&label=Build%20and%20Test
   :alt: GitHub Actions Workflow Status
.. image:: https://img.shields.io/github/actions/workflow/status/esa/polyhedral-gravity-model/.github%2Fworkflows%2Fdocs.yml?logo=gitbook&label=Documentation
   :alt: GitHub Actions Workflow Status

..

.. image:: https://img.shields.io/pypi/v/polyhedral-gravity
   :alt: PyPI - Version
.. image:: https://img.shields.io/badge/platform-linux--64_%7C_win--64_%7C_osx--64_%7C_linux--arm64_%7C_osx--arm64-lightgrey
   :alt: Static Badge
.. image:: https://img.shields.io/pypi/dm/polyhedral-gravity
   :alt: PyPI - Downloads

..

.. image:: https://img.shields.io/conda/v/conda-forge/polyhedral-gravity-model
   :alt: Conda Version
.. image:: https://img.shields.io/conda/pn/conda-forge/polyhedral-gravity-model
   :alt: Conda Platform
.. image:: https://img.shields.io/conda/dn/conda-forge/polyhedral-gravity-model
   :alt: Conda Downloads


Have a look at the **INSTALLATION & QUICK START** section to quickly
get into the use of the polyhedral gravity model.

And for more details, refer to the **Python API** or **C++ API**.

If this implementation proves useful to you, please consider citing the
`accompanying paper <https://doi.org/10.21105/joss.06384>`__
published in the *Journal of Open Source Software*:

.. code-block:: bibtex

   @article{Schuhmacher_Efficient_Polyhedral_Gravity_2024,
      author = {Schuhmacher, Jonas and Blazquez, Emmanuel and Gratl, Fabio and Izzo, Dario and Gómez, Pablo},
      doi = {10.21105/joss.06384},
      journal = {Journal of Open Source Software},
      month = jun,
      number = {98},
      pages = {6384},
      title = {{Efficient Polyhedral Gravity Modeling in Modern C++ and Python}},
      url = {https://joss.theoj.org/papers/10.21105/joss.06384},
      volume = {9},
      year = {2024}
   }

.. toctree::
   :caption: INSTALLATION & QUICK START
   :maxdepth: 2

   quickstart/installation
   quickstart/overview
   quickstart/supported_input
   quickstart/examples_python
   quickstart/examples_cpp

.. toctree::
   :caption: BACKGROUND
   :maxdepth: 2

   background/rational
   background/approach
   background/mesh_integrity
   background/evaluable_vs_evaluate

.. toctree::
   :caption: PYTHON API REFERENCE
   :maxdepth: 2

   api/python

.. toctree::
   :caption: C++ API REFERENCE
   :maxdepth: 2

   api/model
   api/input
   api/output
   api/util

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
