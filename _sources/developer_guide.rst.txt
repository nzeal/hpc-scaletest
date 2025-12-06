===============
Developer Guide
===============

This guide helps developers contribute to and extend HPC-ScaleTest.

Development Setup
=================

Clone and Install
-----------------

.. code-block:: bash

   git clone https://github.com/yourusername/hpc-scaletest.git
   cd hpc-scaletest
   python3 -m venv venv
   source venv/bin/activate
   pip install -e ".[dev]"

Run Tests
---------

.. code-block:: bash

   pytest tests/ -v
   pytest tests/ --cov=. --cov-report=html

Code Quality
------------

.. code-block:: bash

   # Format code
   black .
   
   # Lint
   flake8 .
   
   # Type checking
   mypy .

Architecture Overview
=====================

See :doc:`architecture` for detailed architecture documentation.

Core Concepts
-------------

* **Abstraction**: All backends implement abstract base classes
* **Factory Pattern**: Backends created via factory
* **Registry**: Plugin system for custom backends
* **Configuration**: Dataclasses for type-safe config

Adding a New Scheduler
======================

1. Create Backend Class
-----------------------

.. code-block:: python

   # backends/schedulers/pbs.py
   from core.abstracts import SchedulerBackend
   
   class PBSScheduler(SchedulerBackend):
       def submit_job(self, script_path: str) -> str:
           """Submit PBS job and return job ID."""
           result = subprocess.run(
               ["qsub", script_path],
               capture_output=True, text=True
           )
           return result.stdout.strip()
       
       def get_job_status(self, job_id: str) -> str:
           """Get job status."""
           result = subprocess.run(
               ["qstat", job_id],
               capture_output=True, text=True
           )
           # Parse and return status
           return status
       
       def cancel_job(self, job_id: str) -> bool:
           """Cancel job."""
           result = subprocess.run(
               ["qdel", job_id]
           )
           return result.returncode == 0

2. Register Backend
-------------------

.. code-block:: python

   # In backends/schedulers/__init__.py
   from .pbs import PBSScheduler
   
   __all__ = ['PBSScheduler']

3. Add to Factory
-----------------

.. code-block:: python

   # core/factory.py
   from backends.schedulers.pbs import PBSScheduler
   
   SCHEDULER_MAP = {
       'slurm': SlurmScheduler,
       'pbs': PBSScheduler,  # Add this
       'local': LocalScheduler,
   }

4. Add Tests
------------

.. code-block:: python

   # tests/test_pbs_scheduler.py
   def test_pbs_submit():
       scheduler = PBSScheduler()
       # Test submission logic

Adding a New Launcher
=====================

Similar process to schedulers:

.. code-block:: python

   # backends/launchers/aprun.py
   from core.abstracts import LauncherBackend
   
   class AprunLauncher(LauncherBackend):
       def build_command(self, executable, args, num_tasks, 
                        tasks_per_node=None, **kwargs):
           cmd = ["aprun", "-n", str(num_tasks)]
           if tasks_per_node:
               cmd.extend(["-N", str(tasks_per_node)])
           cmd.append(executable)
           cmd.extend(args)
           return cmd

Custom Input Parsers
====================

For application-specific input formats:

.. code-block:: python

   # utils/custom_parser.py
   from utils.generic_input_parser import GenericInputParser
   
   class FortranNamelistParser(GenericInputParser):
       def parse(self, content: str) -> dict:
           \"\"\"Parse Fortran namelist format.\"\"\"
           params = {}
           # Parse namelist format
           return params
       
       def update(self, content: str, updates: dict) -> str:
           \"\"\"Update namelist parameters.\"\"\"
           # Update and return new content
           return new_content

Testing Guidelines
==================

Unit Tests
----------

.. code-block:: python

   # tests/test_myfeature.py
   import pytest
   from core.test_definition import Test
   
   def test_basic_configuration():
       test = Test(name="test", command=["./app"])
       test.set_resources(max_nodes=4)
       assert test.resource_config.max_nodes == 4

Integration Tests
-----------------

.. code-block:: python

   @pytest.mark.integration
   def test_full_workflow():
       # Test complete workflow
       pass

Mock External Dependencies
--------------------------

.. code-block:: python

   from unittest.mock import patch, MagicMock
   
   @patch('subprocess.run')
   def test_job_submission(mock_run):
       mock_run.return_value = MagicMock(
           stdout="12345", returncode=0
       )
       # Test with mocked subprocess

Documentation Guidelines
========================

Docstring Format
----------------

Use Google-style docstrings:

.. code-block:: python

   def my_function(param1: str, param2: int) -> bool:
       \"\"\"Short description.
       
       Longer description with more details about what
       this function does.
       
       Args:
           param1: Description of param1
           param2: Description of param2
       
       Returns:
           Description of return value
       
       Raises:
           ValueError: When input is invalid
       
       Example:
           >>> my_function("test", 42)
           True
       \"\"\"
       pass

RST Documentation
-----------------

Update relevant .rst files when adding features.

Code Style Guidelines
=====================

* Follow PEP 8
* Use type hints
* Maximum line length: 88 characters (Black default)
* Use descriptive variable names
* Add docstrings to all public functions/classes

Pull Request Process
====================

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Run test suite
6. Update documentation
7. Submit pull request

Release Process
===============

1. Update version in setup.py
2. Update CHANGELOG.md
3. Create Git tag
4. Build and publish to PyPI

.. code-block:: bash

   python setup.py sdist bdist_wheel
   twine upload dist/*

See Also
========

* :doc:`architecture` - Architecture details
* :doc:`contributing` - Contribution guidelines
* :doc:`api_reference` - API documentation
