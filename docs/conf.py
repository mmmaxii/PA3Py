import os
import sys

# Add the source directory to sys.path so autodoc can find the module
sys.path.insert(0, os.path.abspath('../src'))

import pa3py

project = 'PA3Py'
copyright = '2026, Maximiliano Valderrama'
author = 'Maximiliano Valderrama'
release = pa3py.__version__

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'nbsphinx',
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# nbsphinx configuration
nbsphinx_execute = 'never' # Do not execute notebooks during build, use them as they are
