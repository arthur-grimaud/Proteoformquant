# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath("../proteoformquant/Classes"))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Proteoformquant"
copyright = "2023, Arthur Grimaud"
author = "Arthur Grimaud"
release = "0.1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
]

html_theme = "sphinx_rtd_theme"


# mock the folowing modules to avoid errors when building the documentation:
autodoc_mock_imports = [
    "dash",
    "dash_bootstrap_components",
    "dash_daq",
    "fnnls",
    "jsonc_parser",
    "kneed",
    "matplotlib",
    "mpmath",
    "networkx",
    "numpy=",
    "pandas",
    "plotly",
    "progress",
    "pymoo",
    "pyteomics",
    "PyYAML",
    "scikit_learn",
    "scipy",
    "spectrum_utils",
    "SQLAlchemy",
    "sympy",
    "unimod_mapper",
    "lxml",
    "ms_deisotope",
    "numpy",
    "alive_progress",
    "sqlalchemy",
    "Classes",
    "Utils",
    "numba",
    "sklearn",
]
