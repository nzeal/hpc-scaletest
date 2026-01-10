# Configuration file for the Sphinx documentation builder.

import os
import sys
from unittest.mock import MagicMock

# -- Path setup --------------------------------------------------------------
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
project = 'HPC-ScaleTest'
copyright = '2024, HPC-ScaleTest Contributors'
author = 'HPC-ScaleTest Contributors'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Extension configuration -------------------------------------------------

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}

# Suppress warnings
suppress_warnings = ['app.add_node', 'app.add_directive', 'app.add_role']

# -- Mock imports for packages that may not be installed --------------------

class Mock(MagicMock):
    """Enhanced mock for Sphinx autodoc."""
    
    @classmethod
    def __getattr__(cls, name):
        return MagicMock()
    
    def __call__(self, *args, **kwargs):
        return MagicMock()

MOCK_MODULES = [
    'yaml',
    'numpy',
    'matplotlib',
    'matplotlib.pyplot',
    'pandas',
    'scipy',
    'scipy.stats',
    'mpi4py',
    'mpi4py.MPI',
    'backends',
    'backends.schedulers',
    'backends.schedulers.slurm',
    'backends.schedulers.local',
    'backends.launchers',
    'backends.launchers.srun',
    'backends.launchers.mpirun',
    'backends.launchers.simple',
    'backends.modules',
    'backends.modules.lmod',
    'backends.modules.tmod',
    'backends.modules.tmod4',
    'backends.modules.nomod',
    'backends.builds',
    'backends.builds.cmake',
    'backends.builds.make',
    'backends.builds.autotools',
    'backends.builds.spack',
    'backends.builds.easybuild',
    'core',
    'core.test_definition',
    'core.config',
    'core.abstracts',
    'core.factory',
    'core.types',
    'engine',
    'engine.orchestrator',
    'engine.runner',
    'engine.scaling',
    'engine.job_builder',
    'utils',
    'utils.config_parser',
    'utils.generic_input_parser',
    'utils.parsers',
    'utils.code_acquisition',
    'utils.readme_analyzer',
    'utils.input_analyzer',
    'utils.parameter_suggestion',
    'utils.llm_parameter_mapper',
    'utils.input_generator',
    'utils.input_file_scaler',
    'utils.report_generator',
    'utils.scaling_visualizer',
    'utils.system_info',
    'utils.system_loader',
    'utils.slurm_detector',
    'utils.job_submitter',
    'utils.validators',
    'utils.logging_config',
    'utils.file_utils',
]

sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

# -- MathJax configuration ---------------------------------------------------
mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'

mathjax3_config = {
    'tex': {
        'inlineMath': [['$', '$'], ['\\(', '\\)']],
        'displayMath': [['$$', '$$'], ['\\[', '\\]']],
    }
}
