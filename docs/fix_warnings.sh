#!/bin/bash
# Quick fix script for Sphinx build warnings

set -e

echo "================================================"
echo "  HPC-ScaleTest Documentation - Warning Fix"
echo "================================================"
echo ""

# Check if we're in the right directory
if [ ! -f "source/conf.py" ]; then
    echo "❌ Error: Please run this script from the docs_output directory"
    echo "   Usage: cd docs_output && bash fix_warnings.sh"
    exit 1
fi

echo "✓ Found documentation directory"
echo ""

# Backup original conf.py
echo "📋 Backing up original conf.py..."
cp source/conf.py source/conf.py.backup
echo "✓ Backup saved to source/conf.py.backup"
echo ""

# Create updated conf.py with mock imports
echo "🔧 Creating updated conf.py with mock imports..."
cat > source/conf.py << 'EOFCONF'
# Configuration file for the Sphinx documentation builder.

import os
import sys
from pathlib import Path

# Add project root to path for autodoc
sys.path.insert(0, str(Path('..', '..').resolve()))

# Project information
project = 'HPC-ScaleTest'
copyright = '2024, HPC ScaleTest Contributors'
author = 'HPC ScaleTest Contributors'
release = '1.0.0'
version = '1.0'

# General configuration
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.ifconfig',
    'sphinx.ext.githubpages',
]

templates_path = ['_templates']
exclude_patterns = []
source_suffix = '.rst'
master_doc = 'index'

# HTML output
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_theme_options = {
    'canonical_url': '',
    'logo_only': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'style_nav_header_background': '#2980B9',
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# Mock imports for autodoc
autodoc_mock_imports = [
    'core',
    'engine',
    'backends',
    'utils',
    'yaml',
]

# Autodoc options
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
}

# Todo extension
todo_include_todos = True

# Suppress minor warnings
suppress_warnings = ['ref.citation']
EOFCONF

echo "✓ Updated conf.py created"
echo ""

# Clean previous build
echo "🧹 Cleaning previous build..."
make clean > /dev/null 2>&1
echo "✓ Build directory cleaned"
echo ""

# Rebuild documentation
echo "🔨 Building HTML documentation..."
echo ""
if make html; then
    echo ""
    echo "================================================"
    echo "  ✅ Success! Documentation built with minimal warnings"
    echo "================================================"
    echo ""
    echo "📂 HTML documentation: build/html/index.html"
    echo ""
    echo "View locally:"
    echo "  cd build/html && python3 -m http.server 8000"
    echo ""
    echo "To restore original conf.py:"
    echo "  cp source/conf.py.backup source/conf.py"
    echo ""
else
    echo ""
    echo "❌ Build had errors. Restoring backup..."
    cp source/conf.py.backup source/conf.py
    echo "✓ Original conf.py restored"
    exit 1
fi
