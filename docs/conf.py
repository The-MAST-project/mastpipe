"""Sphinx configuration for the mastpipe documentation."""

from __future__ import annotations

from mastpipe import __version__

project = "mastpipe"
author = "The MAST project"
copyright = "2026, The MAST project"
release = __version__
version = __version__

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "sphinx_design",
]

myst_enable_extensions = ["colon_fence", "deflist", "fieldlist", "substitution"]
myst_heading_anchors = 3

# Tutorial notebooks are committed without outputs and executed at build time.
nb_execution_mode = "cache"
nb_execution_timeout = 900
nb_execution_raise_on_error = True
nb_merge_streams = True

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented"
always_document_param_types = False
typehints_fully_qualified = False
napoleon_numpy_docstring = True
napoleon_google_docstring = True
autosummary_generate = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "astropy": ("https://docs.astropy.org/en/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
}
# Offline builds (CI without network) still succeed: intersphinx failures are warnings
# only when the inventory cannot be fetched, and we keep them from failing -W.
suppress_warnings = [
    "intersphinx.external",
    "myst.header",
    "ref.python",
    "sphinx_autodoc_typehints.forward_reference",  # pydantic's own Field annotations
]
nitpicky = False

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "design",
    "jupyter_execute",
    "**.ipynb_checkpoints",
]

html_theme = "furo"
html_title = f"mastpipe {release}"
html_static_path = ["_static"]
html_theme_options = {
    "source_repository": "",
    "navigation_with_keys": True,
}
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
