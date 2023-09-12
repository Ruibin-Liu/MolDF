# Configuration file for the Sphinx documentation builder.
import sys

sys.path.append("..")

# -- Project information

project = "pdbx2df"
copyright = "2023, Ruibin Liu"
author = "Ruibin Liu"

release = "0.6.2"
version = "0.6.2"

# -- General configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]

# -- Options for HTML output

html_theme = "sphinx_rtd_theme"

# -- Options for EPUB output
# epub_show_urls = 'footnote'
