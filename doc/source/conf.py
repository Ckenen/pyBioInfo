# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys

# Point autodoc at the package source tree (repo-root/src), so the
# documentation can be built without installing the package locally.
sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 'src')))


# -- Project information -----------------------------------------------------

project = 'pyBioInfo'
copyright = '2022, Zonggui Chen'
author = 'Zonggui Chen'

# The full version, including alpha/beta/rc tags
release = '1.0.1'


# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # supports NumPy / Google style docstrings
    'sphinx.ext.viewcode',  # links documented objects to source code
    'sphinx.ext.intersphinx',  # links to the Python standard library docs
]

# Autodoc options: keep the docstrings of __init__ methods.
autoclass_content = 'both'

# Intersphinx mapping to the Python standard library.
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
}

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
