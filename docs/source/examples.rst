========
Examples
========

This page provides complete configuration examples for common use cases.

Strong Scaling Example
======================

.. code-block:: yaml

   # strong_scaling.yaml
   repository: /path/to/code
   scaling:
     type: strong
     nodes: 64
     initial_procs: [2, 2, 2]
   hardware:
     procs_per_node: 128
   scheduler: slurm
   partition: standard

Weak Scaling Example (2D)
==========================

.. code-block:: yaml

   # weak_scaling_2d.yaml
   repository: /path/to/code
   scaling:
     type: weak
     nodes: 32
     scaling_factor: 2
     scaling_dimensions: 2
   hardware:
     procs_per_node: 128

GPU Example
===========

.. code-block:: yaml

   # gpu_scaling.yaml
   hardware:
     type: gpu
     gpus_per_node: 4
   scheduler: slurm
   partition: gpu_queue

See :doc:`user_guide` for detailed examples.
