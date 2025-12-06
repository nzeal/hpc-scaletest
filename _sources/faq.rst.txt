===
FAQ
===

Frequently Asked Questions about HPC-ScaleTest.

General Questions
=================

What is HPC-ScaleTest?
----------------------

HPC-ScaleTest is a framework for automating scaling benchmarks on HPC systems.

Who should use it?
------------------

Researchers, HPC application developers, and system administrators who need to evaluate application scalability.

What makes it different?
-------------------------

* Fully automated workflow
* Works with any HPC application
* Intelligent parameter detection
* Multi-backend support

Usage Questions
===============

Can I use it on my laptop?
---------------------------

Yes! Use ``scheduler: local`` for testing.

Do I need to modify my code?
-----------------------------

No. HPC-ScaleTest works with existing applications.

What schedulers are supported?
-------------------------------

* SLURM
* PBS/Torque
* Local execution

Can I test GPU applications?
-----------------------------

Yes. Set ``hardware: type: gpu`` in configuration.

Does it work with Python codes?
--------------------------------

Yes. Specify Python command in configuration.

Configuration Questions
=======================

How do I specify custom modules?
---------------------------------

.. code-block:: yaml

   modules:
     - compiler/gcc/11.2.0
     - mpi/openmpi/4.1.1

Can I use environment variables?
---------------------------------

.. code-block:: yaml

   environment:
     OMP_NUM_THREADS: "4"
     TMPDIR: "/scratch"

How do I scale only in 2D?
---------------------------

.. code-block:: yaml

   scaling:
     scaling_dimensions: 2  # X→Y pattern

Troubleshooting
===============

Jobs fail immediately
---------------------

Check:

1. Module loads: ``module list``
2. Executable path
3. Resource limits
4. Input file format

Results look wrong
------------------

Verify:

1. Timing extraction is correct
2. Problem size actually scaled
3. No I/O bottlenecks

Module not found
----------------

Update module names or contact HPC support.

Technical Questions
===================

Can I extend the framework?
----------------------------

Yes! See :doc:`developer_guide`.

Is it open source?
------------------

Yes, MIT licensed.

Can I contribute?
-----------------

Yes! See :doc:`contributing`.

How do I report bugs?
---------------------

Open an issue on GitHub.

Where can I get help?
---------------------

* Documentation: This site
* GitHub Issues
* Email: support@example.com

Advanced Questions
==================

Can I use custom input formats?
--------------------------------

Yes. Override ``get_input_content()`` method.

Can I add custom metrics?
--------------------------

Yes. Extend ``ReportGenerator`` class.

Can I use a different build system?
------------------------------------

Yes. Implement ``BuildBackend`` interface.

How is weak scaling calculated?
--------------------------------

Domain size scales proportionally with process count to maintain constant work per process.

What is 2D vs 3D scaling?
--------------------------

* 2D: Scales X→Y→X→Y
* 3D: Scales X→Y→Z→X→Y→Z

Can I specify custom scaling patterns?
---------------------------------------

Yes. Modify ``engine/scaling.py``.

See Also
========

* :doc:`troubleshooting` - Troubleshooting guide
* :doc:`user_guide` - User documentation
* :doc:`examples` - Example configurations
