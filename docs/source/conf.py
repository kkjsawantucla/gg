import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join('..', '..')))

project = 'graph_gcbh'
copyright = '2024, Kaustubh Sawant'
author = 'Kaustubh Sawant'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
extensions = [
    # other extensions
    'sphinx.ext.githubpages',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx.ext.autodoc',
]

todo_include_todos = True
autoclass_content = 'init'

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
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
