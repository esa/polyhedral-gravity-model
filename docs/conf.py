# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

# -- Read-the-docs setup -----------------------------------------------------
# Adapted from https://github.com/TartanLlama/cpp-documentation-example/blob/master/docs/conf.py
import subprocess, os


def configure_doxyfile(input_dir, output_dir):
    with open("Doxyfile.in", "r") as file:
        filedata = file.read()

    filedata = filedata.replace("@DOXYGEN_INPUT_DIR@", input_dir)
    filedata = filedata.replace("@DOXYGEN_OUTPUT_DIR@", output_dir)

    with open("Doxyfile", "w") as file:
        file.write(filedata)


# If you build the docs locally, but not via CMake, first set the environment variable BUILD_DOCS_CLI
breathe_projects = {}
if any([
    os.environ.get("READTHEDOCS", None),
    os.environ.get("GITHUB_PAGES_BUILD", None),
    os.environ.get("BUILD_DOCS_CLI", None)]
):
    # These lines assume that we are in  /<repo-root>/docs/ folder
    input_dir = "../src"
    output_dir = "build"
    configure_doxyfile(input_dir, output_dir)
    subprocess.call("doxygen", shell=True)
    breathe_projects["polyhedral-gravity-model"] = output_dir + "/xml"

# -- Project information -----------------------------------------------------

project = "Polyhedral Gravity Model"
copyright = "2022 - 2025, Jonas Schuhmacher"
author = "Jonas Schuhmacher"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["breathe", "sphinx.ext.napoleon", "sphinx.ext.autodoc"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static", "figures/eros_010.png", "figures/runtime_measurements.png"]

# Breathe Configuration
breathe_default_project = "polyhedral-gravity-model"
breathe_default_members = ("members", "undoc-members")

# -- Options for Theme 'sphinx_book_theme' -----------------------------------
# https://sphinx-book-theme.readthedocs.io/en/latest/tutorials/get-started.html
html_theme_options = {
    "use_edit_page_button": True,
    "use_source_button": True,
    "use_issues_button": True,
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/esa/polyhedral-gravity-model",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        },
        {
            "name": "Conda Forge",
            "url": "https://anaconda.org/conda-forge/polyhedral-gravity-model",
            "icon": "https://img.shields.io/conda/vn/conda-forge/polyhedral-gravity-model",
            "type": "url",
        },
        {
            "name": "PyPi",
            "url": "https://pypi.org/project/polyhedral-gravity/",
            "icon": "https://img.shields.io/pypi/v/polyhedral-gravity",
            "type": "url",
        },
    ],
}

html_context = {
    "github_url": "https://github.com",
    "github_user": "esa",
    "github_repo": "polyhedral-gravity-model",
    "github_version": "main",
    "doc_path": "docs",
}
