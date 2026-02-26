from datetime import datetime
from pathlib import Path
import sys

project = 'PyMKMKit'
author = 'PyMKMKit contributors'
copyright = f"{datetime.now():%Y}, {author}"

# Ensure autodoc can import the local package when docs are built from the
# repository checkout.
ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

autodoc_member_order = 'bysource'
autodoc_typehints = 'description'
autosummary_generate = True

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/logo/logo-256.png'
html_favicon = '_static/logo/favicon.ico'
html_css_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.1.1/css/all.min.css"
]