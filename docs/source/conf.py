# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'graph_gcbh'
copyright = '2024, Kaustubh Sawant'
author = 'Kaustubh Sawant'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configurations

extensions = [
    # other extensions
    'sphinx.ext.githubpages',
    'sphinx.ext.autodoc',
]

html_context = {
    "display_github": True,  # Integrate GitHub
    "github_user": "kkjsawantucla",  # Username
    "github_repo": "gg",  # Repo name
    "github_version": "main",  # Version
    "conf_py_path": "/docs/source/",  # Path in the checkout to the docs root
}

templates_path = ['_templates']
exclude_patterns = []

language = 'English'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
