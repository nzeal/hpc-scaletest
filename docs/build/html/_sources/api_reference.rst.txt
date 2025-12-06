=============
API Reference
=============

Complete API documentation for HPC-ScaleTest.

Core Modules
============

.. toctree::
   :maxdepth: 2

   api/core
   api/engine
   api/backends
   api/utils

Quick Reference
===============

Test Definition
---------------

.. code-block:: python

   from core.test_definition import Test
   test = Test(name="my_test", command=["./app"])

Orchestrator
------------

.. code-block:: python

   from engine.orchestrator import HPCOrchestrator
   orchestrator = HPCOrchestrator(config)
   orchestrator.run()

See module documentation for complete API details.
