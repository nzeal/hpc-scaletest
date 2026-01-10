============
Architecture
============

This document describes the internal architecture of HPC-ScaleTest.

System Overview
===============

HPC-ScaleTest follows a layered architecture:

.. code-block:: text

   ┌─────────────────────────────────────┐
   │     User Interface Layer            │
   │  (CLI, Python API, YAML Config)     │
   └─────────────────────────────────────┘
                    │
   ┌─────────────────────────────────────┐
   │     Orchestration Layer             │
   │  (Orchestrator, TestRunner)         │
   └─────────────────────────────────────┘
                    │
   ┌─────────────────────────────────────┐
   │     Core Layer                      │
   │  (Test Definition, Config, Types)   │
   └─────────────────────────────────────┘
                    │
   ┌─────────────────────────────────────┐
   │     Backend Abstraction Layer       │
   │  (Abstract Base Classes, Factory)   │
   └─────────────────────────────────────┘
                    │
   ┌─────────────────────────────────────┐
   │     Backend Implementations         │
   │  (Schedulers, Launchers, Modules)   │
   └─────────────────────────────────────┘

Core Components
===============

Test Definition (core/)
-----------------------

Provides user-facing API for defining tests.

**Key Classes**:

* ``Test``: Main test configuration class
* ``BackendConfig``: Backend selection
* ``ResourceConfig``: Resource allocation
* ``ScalingConfig``: Scaling parameters

Orchestration (engine/)
------------------------

Manages workflow execution.

**Key Classes**:

* ``HPCOrchestrator``: Complete workflow automation
* ``TestRunner``: Job execution and monitoring
* ``ScalingEngine``: Scaling algorithm implementation
* ``JobBuilder``: Job script generation

Backends (backends/)
--------------------

Pluggable implementations for different systems.

**Schedulers**: SLURM, PBS, Local
**Launchers**: srun, mpirun, aprun
**Modules**: Lmod, Tmod, no-modules
**Build Systems**: CMake, Make, Autotools, Spack

Utilities (utils/)
------------------

Helper functions and tools.

**Categories**:

* Configuration parsing
* Code analysis
* Input file generation
* Report generation
* System detection

Design Patterns
===============

Abstract Factory Pattern
------------------------

Used for creating backend instances:

.. code-block:: python

   # Abstract products
   class SchedulerBackend(ABC):
       @abstractmethod
       def submit_job(self, script): pass
   
   # Concrete products
   class SlurmScheduler(SchedulerBackend):
       def submit_job(self, script):
           # SLURM implementation
   
   # Factory
   class BackendFactory:
       def create_scheduler(self, name):
           return SCHEDULER_MAP[name]()

Registry Pattern
----------------

For plugin system:

.. code-block:: python

   class Registry:
       _backends = {}
       
       @classmethod
       def register(cls, name, backend_class):
           cls._backends[name] = backend_class
       
       @classmethod
       def get(cls, name):
           return cls._backends[name]

Strategy Pattern
----------------

For scaling algorithms:

.. code-block:: python

   class ScalingStrategy(ABC):
       @abstractmethod
       def generate_configs(self, base_config):
           pass
   
   class StrongScaling(ScalingStrategy):
       def generate_configs(self, base_config):
           # Strong scaling logic
   
   class WeakScaling(ScalingStrategy):
       def generate_configs(self, base_config):
           # Weak scaling logic

Data Flow
=========

Configuration Flow
------------------

.. code-block:: text

   YAML File → ConfigParser → Config Objects → Test
                    ↓
   Python API → Direct Config → Config Objects → Test
                    ↓
   CLI Args → ArgParser → Config Objects → Test

Execution Flow
--------------

.. code-block:: text

   1. Code Acquisition
      ↓
   2. README Analysis → Build Info
      ↓
   3. Code Compilation → Executable
      ↓
   4. Input Analysis → Parameters
      ↓
   5. Scaling Engine → Job Configs
      ↓
   6. Job Builder → Job Scripts
      ↓
   7. Job Submission → Job IDs
      ↓
   8. Job Monitoring → Status Updates
      ↓
   9. Results Collection → Raw Data
      ↓
   10. Report Generation → Metrics

Module Dependencies
===================

.. code-block:: text

   core.test_definition
       ├── core.config
       ├── core.types
       └── core.abstracts
   
   engine.orchestrator
       ├── core.test_definition
       ├── engine.runner
       ├── utils.code_acquisition
       ├── utils.readme_analyzer
       └── utils.report_generator
   
   engine.runner
       ├── core.test_definition
       ├── engine.scaling
       ├── engine.job_builder
       └── backends.*
   
   backends.*
       └── core.abstracts

Extension Points
================

Custom Schedulers
-----------------

Implement ``SchedulerBackend`` abstract class.

Custom Launchers
----------------

Implement ``LauncherBackend`` abstract class.

Custom Input Parsers
--------------------

Inherit from ``GenericInputParser``.

Custom Scaling Algorithms
--------------------------

Implement scaling logic in ``engine/scaling.py``.

Error Handling
==============

Exception Hierarchy
-------------------

.. code-block:: text

   HPCScaleTestError (base)
       ├── ConfigurationError
       ├── ValidationError
       ├── BuildError
       ├── JobSubmissionError
       └── ResultsCollectionError

Error Propagation
-----------------

Errors are caught at the orchestrator level and logged appropriately.

Performance Considerations
==========================

* Lazy loading of modules
* Caching of system detection
* Parallel job submission (future)
* Efficient file I/O

Security Considerations
=======================

* No hardcoded credentials
* Safe subprocess execution
* Input sanitization
* File permission checks

Testing Strategy
================

Unit Tests
----------

Test individual components in isolation.

Integration Tests
-----------------

Test component interactions.

System Tests
------------

End-to-end workflow tests.

See Also
========

* :doc:`developer_guide` - Development guide
* :doc:`api_reference` - API documentation
