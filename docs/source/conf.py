# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os 
import sys
import click

# -- Path setup --------------------------------------------------------------
# Since the package lives under src/, we need to add it to sys.path so that
# autodoc can import the modules and extract docstrings.
sys.path.insert(0, os.path.abspath('../../src'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'oddSNP'
copyright = '2026, Rodolfo Allendes'
author = 'Rodolfo Allendes'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
  'sphinx.ext.autodoc',  # auto-generate from docstrings
  'sphinx.ext.napoleon', # Google/NumPy style docstrings
  'sphinx.ext.viewcode', # add links to source code
  'sphinx.ext.autosummary', # generate summary tables for modules/classes/functions
  'nbsphinx', # include Jupyter notebooks in the documentation
  'sphinx_click', # support for documenting Click CLI commands
]

# Napoleon settings to support Google style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True

# Autodoc settings
autodoc_default_options = {
  'members': True, # include members of modules/classes
  'undoc-members': False, # include members without docstrings
  'private-members': False, # include private members (starting with _)
  'show-inheritance': True, # show class inheritance diagrams
  'member-order': 'alphabetical', # order members as they appear in source code
}
autodoc_typehints = 'description' # show type hints in the description instead of signature
autosummary_generate = True # generate stub files for autosummary

# nbsphinx settings
nbsphinx_execute = 'never' # do not execute notebooks when building docs
nbsphinx_allow_errors = True # allow errors in notebooks without breaking the build

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
html_css_files = ['custom.css']
html_title = 'oddSNP'

html_logo = '_static/logo.png'

html_theme_options = {
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
}
