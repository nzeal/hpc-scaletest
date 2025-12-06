=============
Configuration
=============

Configuration guide for HPC-ScaleTest.

YAML Configuration
==================

Basic configuration file structure:

.. code-block:: yaml

   repository: /path/to/code
   scaling:
     type: weak
     nodes: 32
   hardware:
     procs_per_node: 128
   scheduler: slurm
   partition: standard

Python API Configuration
========================

.. code-block:: python

   test = Test(name="my_test")
   test.set_backend(scheduler="slurm")
   test.set_resources(max_nodes=32)
   test.set_scaling(scaling_type="weak")

See :doc:`user_guide` for detailed configuration options.
