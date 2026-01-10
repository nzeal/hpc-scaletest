===============
Troubleshooting
===============

Common Issues and Solutions
===========================

Installation Issues
===================

ImportError: No module named 'yaml'
------------------------------------

**Solution**: Install PyYAML

.. code-block:: bash

   pip install pyyaml

Permission Denied
-----------------

**Solution**: Use --user flag

.. code-block:: bash

   pip install --user -r requirements.txt

Job Submission Issues
=====================

Jobs Fail Immediately
---------------------

**Symptoms**: Jobs exit right after submission

**Diagnosis**:

1. Check SLURM output: ``cat slurm-XXXXX.out``
2. Check error log: ``cat slurm-XXXXX.err``

**Common Causes**:

* Module load failures
* Incorrect executable path
* Missing dependencies
* Insufficient resources

**Solutions**:

.. code-block:: bash

   # Test module loads
   module purge
   module load gcc/11.2.0 openmpi/4.1.1
   
   # Test executable
   ./app --help
   
   # Check resource limits
   ulimit -a

Module Not Found
----------------

**Solution**: Update module names

.. code-block:: bash

   module avail  # List available modules
   module spider gcc  # Search for modules

SLURM: Invalid Account
-----------------------

**Solution**: Check valid accounts

.. code-block:: bash

   sacctmgr show assoc where user=$USER

Results Issues
==============

Missing Timing Data
-------------------

**Symptoms**: Report shows "N/A" for times

**Solutions**:

1. Verify output files exist
2. Check timing format in output
3. Update timing parser if needed

Incorrect Efficiency Metrics
-----------------------------

**Symptoms**: Efficiency > 100% or negative

**Causes**:

* Timing variation
* Wrong baseline
* Incorrect parsing

**Solutions**:

1. Run multiple trials
2. Verify baseline configuration
3. Check timing extraction logic

Build Issues
============

CMake Configuration Failed
---------------------------

**Solution**: Specify compilers

.. code-block:: yaml

   build:
     flags:
       CMAKE_CXX_COMPILER: mpicxx
       CMAKE_C_COMPILER: mpicc

Missing Dependencies
--------------------

**Solution**: Load required modules

.. code-block:: yaml

   modules:
     - compiler/gcc/11.2.0
     - libraries/hdf5/1.12.0

Configuration Issues
====================

YAML Parsing Error
------------------

**Symptoms**: "YAML parsing failed"

**Solution**: Check YAML syntax

.. code-block:: bash

   python -c "import yaml; yaml.safe_load(open('run.yaml'))"

Invalid Configuration
---------------------

**Symptoms**: "Configuration validation failed"

**Solution**: Run validator

.. code-block:: bash

   python validate_scaling.py --config run.yaml

Performance Issues
==================

Poor Scaling Efficiency
-----------------------

**Possible Causes**:

* Communication overhead
* Load imbalance
* I/O contention
* Memory bandwidth

**Diagnosis Tools**:

.. code-block:: bash

   # Profile application
   module load vtune
   srun --profile=vtune ./app

Slow Job Submission
-------------------

**Solution**: Submit in batches

.. code-block:: yaml

   # Limit concurrent submissions
   max_concurrent_jobs: 10

System-Specific Issues
======================

SLURM: Job Pending Forever
---------------------------

**Check**: Queue status and limits

.. code-block:: bash

   squeue -u $USER
   scontrol show partition
   sacctmgr show qos

PBS: Job Not Starting
----------------------

**Check**: PBS status

.. code-block:: bash

   qstat -f $JOBID
   qstat -Q

Local: Out of Memory
--------------------

**Solution**: Reduce problem size or process count

.. code-block:: yaml

   hardware:
     procs_per_node: 4  # Reduce

Getting More Help
=================

Enable Debug Logging
--------------------

.. code-block:: python

   import logging
   logging.basicConfig(level=logging.DEBUG)

Collect Diagnostic Information
-------------------------------

.. code-block:: bash

   # System info
   uname -a
   python --version
   module list
   
   # Configuration
   cat run.yaml
   
   # Error logs
   cat slurm-*.err
   cat output/logs/error.log

Report Issues
-------------

When reporting issues, include:

1. HPC-ScaleTest version
2. System information
3. Configuration files
4. Error messages
5. Steps to reproduce

Contact Support
---------------

* GitHub Issues: https://github.com/user/hpc-scaletest/issues
* Email: support@example.com
* Documentation: https://hpc-scaletest.readthedocs.io

See Also
========

* :doc:`faq` - Frequently asked questions
* :doc:`user_guide` - User documentation
* :doc:`configuration` - Configuration options
