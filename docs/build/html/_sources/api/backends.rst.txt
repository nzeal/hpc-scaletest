===============
Backends Module
===============

The backends module provides interfaces to various HPC system components.

Schedulers
==========

Job scheduling system interfaces.

**Supported Schedulers:**

- **SLURM**: Most common HPC scheduler
- **PBS/Torque**: Traditional HPC scheduler
- **Local**: For local testing without a scheduler

**Usage Example:**

.. code-block:: python

   from backends.schedulers import SlurmScheduler
   
   scheduler = SlurmScheduler()
   job_id = scheduler.submit_job("job_script.sh")

Launchers
=========

MPI job launcher interfaces.

**Supported Launchers:**

- **srun**: SLURM launcher
- **mpirun**: OpenMPI/MPICH launcher
- **mpiexec**: Standard MPI launcher
- **simple**: Direct execution (for testing)

**Usage Example:**

.. code-block:: python

   from backends.launchers import SrunLauncher
   
   launcher = SrunLauncher()
   launcher.run(num_procs=64, executable="./myapp")

Module System
=============

Environment module system interfaces.

**Supported Systems:**

- **Lmod**: Modern Lua-based module system
- **Tmod**: Traditional Tcl module system
- **Tmod4**: Environment Modules 4.x
- **nomod**: No module system (manual environment)

**Usage Example:**

.. code-block:: python

   from backends.modules import LmodModules
   
   modules = LmodModules()
   modules.load("gcc/11.2.0")
   modules.load("openmpi/4.1.0")

Build Systems
=============

Application build system interfaces.

**Supported Systems:**

- **CMake**: Modern build system
- **Make**: Traditional make-based builds
- **Autotools**: GNU autotools (configure/make)
- **Spack**: HPC package manager
- **EasyBuild**: HPC software build framework

**Usage Example:**

.. code-block:: python

   from backends.builds import CMakeBuilder
   
   builder = CMakeBuilder()
   builder.configure(source_dir="/path/to/source")
   builder.build()
   builder.install()
