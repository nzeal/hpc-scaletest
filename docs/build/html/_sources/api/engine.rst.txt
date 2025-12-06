=============
Engine Module
=============

The engine module contains the main orchestration and execution logic.

Module: engine.orchestrator
============================

Orchestrates the entire scaling study workflow.

**Key Classes:**

- ``HPCOrchestrator``: Main orchestrator class
- ``OrchestratorConfig``: Configuration for orchestrator

**Usage Example:**

.. code-block:: python

   from engine.orchestrator import HPCOrchestrator
   
   orchestrator = HPCOrchestrator(config)
   orchestrator.run_study()

Module: engine.runner
=====================

Executes individual test runs.

**Key Classes:**

- ``TestRunner``: Runs individual tests
- ``RunConfig``: Configuration for test runs

Module: engine.scaling
======================

Implements scaling algorithms (weak and strong scaling).

**Key Classes:**

- ``WeakScaling``: Weak scaling configuration generator
- ``StrongScaling``: Strong scaling configuration generator
- ``MPITopologyCalculator``: Calculates optimal MPI topologies

**See detailed documentation in:**
- :doc:`../algorithms/weak_scaling_complete`
- :doc:`../algorithms/strong_scaling_complete`
- :doc:`../algorithms/mpi_topology_complete`

Module: engine.job_builder
===========================

Builds job submission scripts.

**Key Classes:**

- ``JobBuilder``: Creates job scripts for different schedulers
- ``SlurmJobBuilder``: SLURM-specific job script builder
- ``PBSJobBuilder``: PBS-specific job script builder
