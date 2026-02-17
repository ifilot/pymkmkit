from datetime import datetime

project = 'PyMKMKit'
author = 'PyMKMKit contributors'
copyright = f"{datetime.now():%Y}, {author}"

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/logo/logo-256.png'
html_favicon = '_static/logo/favicon.ico'
