=============
CLI Reference
=============

Command-Line Interface
======================

HPC-ScaleTest provides a command-line interface for quick access to common workflows.

Main Command: hpc_auto.py
==========================

Basic Syntax
------------

.. code-block:: bash

   python hpc_auto.py [SOURCE] [OPTIONS]
   python hpc_auto.py --config CONFIG_FILE

Arguments
---------

Positional Arguments
^^^^^^^^^^^^^^^^^^^^

``SOURCE``
   Path to source code directory or Git repository URL

Options
-------

Scaling Options
^^^^^^^^^^^^^^^

``--scaling {strong,weak}``
   Scaling type (default: strong)

``--nodes N``
   Maximum number of nodes (default: 4)

``--scaling-factor F``
   Scaling factor for weak scaling (default: 2)

``--scaling-dimensions {1,2,3}``
   Scaling dimensions: 1=X-only, 2=X↔Y, 3=X↔Y↔Z (default: 2)

``--initial-procs X,Y,Z``
   Initial process decomposition (default: 2,2,2)

Hardware Options
^^^^^^^^^^^^^^^^

``--hardware {cpu,gpu}``
   Hardware type (default: cpu)

``--procs-per-node N``
   Processes per node (default: 128)

``--gpus-per-node N``
   GPUs per node (default: 0)

Scheduler Options
^^^^^^^^^^^^^^^^^

``--scheduler {slurm,pbs,local}``
   Job scheduler (default: slurm)

``--launcher {srun,mpirun,simple}``
   MPI launcher (default: srun)

``--partition NAME``
   Partition/queue name

``--account NAME``
   Account/project name

``--time-limit HH:MM:SS``
   Time limit (default: 02:00:00)

Configuration Options
^^^^^^^^^^^^^^^^^^^^^

``--config FILE``
   Load configuration from YAML file

``--output-dir DIR``
   Output directory (default: ./output)

``--auto-submit``
   Automatically submit jobs (default: True)

``--no-auto-submit``
   Generate scripts without submitting

Build Options
^^^^^^^^^^^^^

``--build-system {cmake,make,autotools,spack}``
   Build system to use

``--modules MODULE [MODULE ...]``
   Environment modules to load

Examples
========

Basic Strong Scaling
--------------------

.. code-block:: bash

   python hpc_auto.py /path/to/code \\
       --scaling strong \\
       --nodes 32 \\
       --partition standard \\
       --account proj123

Weak Scaling from Git
---------------------

.. code-block:: bash

   python hpc_auto.py https://github.com/user/app.git \\
       --scaling weak \\
       --nodes 64 \\
       --scaling-factor 2 \\
       --scaling-dimensions 2

GPU Scaling
-----------

.. code-block:: bash

   python hpc_auto.py /path/to/gpu-code \\
       --hardware gpu \\
       --gpus-per-node 4 \\
       --partition gpu_queue

Using Configuration File
------------------------

.. code-block:: bash

   python hpc_auto.py --config my_test.yaml

Report Generator
================

Generate reports from existing results:

.. code-block:: bash

   python -m utils.report_generator \\
       --input output/results \\
       --output efficiency_report.txt

Validation Tool
===============

Validate configuration before running:

.. code-block:: bash

   python validate_scaling.py --config run.yaml

Exit Codes
==========

* ``0`` - Success
* ``1`` - Configuration error
* ``2`` - Build failure
* ``3`` - Job submission failure
* ``4`` - Results collection failure

See Also
========

* :doc:`user_guide` - Detailed usage guide
* :doc:`configuration` - Configuration options
* :doc:`examples` - Example configurations
