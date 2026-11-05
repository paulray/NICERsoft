# Configuration file for the Sphinx documentation builder.
#
# For the full list of options see
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

import matplotlib

matplotlib.use("Agg")

# Make the nicer package importable without installing it, and let scripts
# that import sibling modules from scripts/ (e.g. InteractiveLC) find them.
DOCS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(DOCS_DIR)
sys.path.insert(0, REPO_DIR)
sys.path.append(os.path.join(REPO_DIR, "scripts"))

# -- Project information -----------------------------------------------------

project = "NICERsoft"
copyright = "2017-2026, NICERsoft contributors"
author = "NICERsoft contributors"

# -- General configuration ---------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinxarg.ext",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
# Convert quotes and ellipses, but not "--" (which appears in option names)
smartquotes_action = "qe"

myst_enable_extensions = ["colon_fence", "deflist"]
myst_heading_anchors = 3

# These packages cannot be installed on the docs builder (HEASoft/PyXspec) or
# are only needed by a few scripts, so stand-in modules are used instead.
autodoc_mock_imports = [
    "xspec",
    "heasoftpy",
    "cartopy",
    "pyorbital",
    "astroquery",
    "bs4",
    "fuzzywuzzy",
    "pandas",
]

# sphinx-argparse executes each script file (with :filename:) outside of the
# autodoc mock context, so install the mocks for the whole build.
from sphinx.ext.autodoc.mock import MockFinder  # noqa: E402

sys.meta_path.insert(0, MockFinder(autodoc_mock_imports))

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
autodoc_member_order = "bysource"

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_rtype = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "pint": ("https://nanograv-pint.readthedocs.io/en/latest/", None),
}

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_title = "NICERsoft"
html_theme_options = {
    "navigation_depth": 3,
}
html_context = {
    "display_github": True,
    "github_user": "paulray",
    "github_repo": "NICERsoft",
    "github_version": "master",
    "conf_py_path": "/docs/",
}


# -- Check that every script is documented -----------------------------------


def _tracked_scripts():
    """Python scripts in scripts/ (only those tracked by git, if possible)."""
    import subprocess

    try:
        out = subprocess.run(
            ["git", "ls-files", "scripts/*.py"],
            cwd=REPO_DIR,
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        files = [os.path.join(REPO_DIR, f) for f in out.split()]
    except (OSError, subprocess.CalledProcessError):
        files = []
    if not files:
        import glob

        files = glob.glob(os.path.join(REPO_DIR, "scripts", "*.py"))
    return sorted(files)


def check_scripts_documented(app):
    from sphinx.util import logging

    logger = logging.getLogger(__name__)
    pages = os.path.join(app.srcdir, "scripts")
    documented = set()
    for name in os.listdir(pages):
        if name.endswith(".md"):
            with open(os.path.join(pages, name)) as f:
                for line in f:
                    if line.startswith(":filename: ../scripts/"):
                        documented.add(line.split("/")[-1].strip())
    for path in _tracked_scripts():
        with open(path) as f:
            if "def get_parser(" not in f.read():
                continue
        if os.path.basename(path) not in documented:
            logger.warning(
                f"script {os.path.basename(path)} is not documented in docs/scripts/"
            )


def setup(app):
    app.connect("builder-inited", check_scripts_documented)
