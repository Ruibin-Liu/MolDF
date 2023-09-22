# Configuration file for the Sphinx documentation builder.
import os
import sys

sys.path.insert(0, os.path.abspath("../"))

from moldf.version import __version__ as moldf_version  # type: ignore

# -- Project information

project = "moldf"
copyright = "2023, Ruibin Liu"
author = "Ruibin Liu"

version = moldf_version
release = moldf_version

# -- General configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
]

napoleon_google_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]

autodoc_member_order = "bysource"

# -- Options for HTML output

html_theme = "furo"
# html_theme = "sphinx_rtd_theme"

# -- Options for EPUB output
# epub_show_urls = 'footnote'
