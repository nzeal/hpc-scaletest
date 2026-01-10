===========
Core Module
===========

The core module provides fundamental abstractions and APIs for HPC-ScaleTest.

Module: core.test_definition
=============================

Defines test configurations and parameters for HPC scaling studies.

**Key Classes:**

- ``TestDefinition``: Main class for defining test configurations
- ``ScalingConfig``: Configuration for scaling parameters  
- ``ResourceSpec``: Specification of computational resources

**Usage Example:**

.. code-block:: python

   from core.test_definition import TestDefinition
   
   test = TestDefinition(
       name="my_test",
       application="lammps",
       scaling_type="weak"
   )

Module: core.config
===================

Handles configuration file parsing and validation.

**Key Classes:**

- ``Config``: Main configuration class
- ``ConfigValidator``: Validates configuration parameters
- ``ConfigParser``: Parses YAML configuration files

**Example Configuration:**

.. code-block:: yaml

   application:
     name: "lammps"
     version: "2023.08.02"
   
   scaling:
     type: "weak"
     baseline_nodes: 1
     target_nodes: [2, 4, 8]

Module: core.abstracts
======================

Provides abstract base classes for extensibility.

**Key Classes:**

- ``AbstractScheduler``: Base class for job schedulers
- ``AbstractLauncher``: Base class for MPI launchers
- ``AbstractBuilder``: Base class for build systems

Module: core.factory
====================

Implements factory pattern for creating instances.

**Key Functions:**

- ``create_scheduler(name)``: Creates scheduler instance
- ``create_launcher(name)``: Creates launcher instance
- ``create_builder(name)``: Creates builder instance

Module: core.types
==================

Defines type hints and custom types used throughout HPC-ScaleTest.

**Key Types:**

- ``NodeCount``: Type for node counts
- ``ProcessCount``: Type for process counts  
- ``MPITopology``: Type for MPI topology tuples
