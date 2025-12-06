============
Utils Module
============

The utils module provides utility functions and helper classes.

Configuration Utilities
=======================

**Modules:**

- ``utils.config_parser``: Parse configuration files
- ``utils.generic_input_parser``: Generic input file parsing
- ``utils.parsers``: Various parser utilities

Input File Handling
===================

**Modules:**

- ``utils.input_generator``: Generate input files for applications
- ``utils.input_file_scaler``: Scale input files for different process counts
- ``utils.input_analyzer``: Analyze input file structure

**Usage Example:**

.. code-block:: python

   from utils.input_file_scaler import scale_input_file
   
   scale_input_file(
       input_file="input.in",
       scale_factor=2.0,
       output_file="input_scaled.in"
   )

AI-Assisted Utilities
=====================

**Modules:**

- ``utils.readme_analyzer``: Analyze README files for build instructions
- ``utils.parameter_suggestion``: Suggest optimal parameters
- ``utils.llm_parameter_mapper``: Use LLM for parameter mapping

Code Acquisition
================

**Module:** ``utils.code_acquisition``

Download and manage application source code.

**Functions:**

- ``download_from_git(url, dest)``: Clone git repository
- ``download_from_url(url, dest)``: Download tarball
- ``verify_checksum(file, expected)``: Verify file integrity

Reporting and Visualization
============================

**Modules:**

- ``utils.report_generator``: Generate study reports
- ``utils.scaling_visualizer``: Create scaling plots

**Usage Example:**

.. code-block:: python

   from utils.report_generator import generate_report
   
   generate_report(
       results_dir="results/",
       output_file="report.pdf"
   )

System Information
==================

**Modules:**

- ``utils.system_info``: Gather system information
- ``utils.system_loader``: Load system configurations
- ``utils.slurm_detector``: Detect SLURM environment

Job Management
==============

**Module:** ``utils.job_submitter``

Submit and monitor jobs.

**Functions:**

- ``submit_job(script_path)``: Submit a job
- ``check_job_status(job_id)``: Check job status
- ``wait_for_job(job_id)``: Wait for job completion

Validation
==========

**Module:** ``utils.validators``

Validate configurations and inputs.

**Functions:**

- ``validate_config(config)``: Validate configuration
- ``validate_topology(topology, grid)``: Validate MPI topology
- ``validate_paths(paths)``: Validate file paths

Logging
=======

**Module:** ``utils.logging_config``

Configure logging for HPC-ScaleTest.

File Utilities
==============

**Module:** ``utils.file_utils``

File operation utilities.

**Functions:**

- ``copy_file(src, dst)``: Copy files
- ``create_directory(path)``: Create directories
- ``find_files(pattern, directory)``: Find files matching pattern
