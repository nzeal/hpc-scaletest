==========
User Guide
==========

Comprehensive guide to using HPC-ScaleTest.

Concepts
========

Strong Scaling
--------------

Fixed problem size, increasing parallelism.

Weak Scaling
------------

Problem size grows with resources.

Basic Usage
===========

CLI
---

.. code-block:: bash

   python hpc_auto.py /path/to/code --scaling weak --nodes 32

Python API
----------

.. code-block:: python

   from core.test_definition import Test
   test = Test(name="my_test", command=["./app"])
   test.set_scaling(scaling_type="weak", max_nodes=32)

See :doc:`examples` for more usage patterns.
