# HPC-ScaleTest

A modular Python framework for running benchmark scaling tests on heterogeneous HPC systems with support for CPU and GPU resources.

## 🚀 NEW: Automated End-to-End Workflow

**Run complete HPC scaling tests with a single command!**

```bash
# Simple strong scaling test
python hpc_auto.py /path/to/code --scaling strong --nodes 8

# Clone from Git and test
python hpc_auto.py https://github.com/user/hpc-app.git --scaling weak --nodes 16
```

**Features:**
- ✅ **Minimal Input**: Just provide code path or Git URL
- ✅ **Intelligent Analysis**: Auto-detect dependencies from README
- ✅ **Auto-Compilation**: Handles CMake, Make, Autotools automatically
- ✅ **Scaling Tests**: Automated strong/weak scaling configurations
- ✅ **CPU & GPU Support**: Seamless heterogeneous workloads
- ✅ **Efficiency Reports**: Publication-ready performance analysis

➡️ **[See Full Automated Workflow Documentation](docs/AUTOMATED_WORKFLOW.md)**

---

## Core Features

- **Layered Architecture**: Clean separation between user-facing API and backend implementations
- **Scheduler Agnostic**: Write tests once, run on any system (Local, Slurm, PBS, etc.)
- **Pluggable Backends**: Easily swap schedulers, launchers, module systems, and build tools
- **System Configuration**: Define custom launchers and HPC system configurations (ReFrame-style)
- **Scaling Tests**: Built-in support for strong and weak scaling studies
- **Heterogeneous Support**: Handle mixed CPU/GPU resources
- **Automatic Job Generation**: Generate and submit jobs with proper resource allocations
- **Flexible Job Submission**: Control whether jobs are submitted automatically or manually
- **Result Aggregation**: Collect and analyze performance metrics
- **Intelligent Automation**: AI-driven dependency detection and build configuration

## Architecture

```
User Test Definition (Scheduler-Agnostic)
           ↓
    Abstraction Layer (ABCs)
           ↓
Backend Implementations (Pluggable)
  - Schedulers: Local, Slurm, PBS
  - Launchers: srun, mpirun, mpiexec
  - Modules: Lmod, Tmod, None
  - Builds: Make, CMake, Autotools, EasyBuild, Spack
```

## Installation

```bash
# Clone repository
git clone <repo-url>
cd hpc-scaletest

# Install dependencies (minimal requirements)
pip install pyyaml  # Optional, for YAML configs

# Make CLI executable
chmod +x scaletest.py
```

## Quick Start

### Option 1: Automated Workflow (Recommended)

**Simplest way to run scaling tests:**

```bash
# Command-line interface
python hpc_auto.py /path/to/code --scaling strong --nodes 8

# Python API
from engine.orchestrator import create_simple_workflow

orchestrator = create_simple_workflow(
    source="/path/to/code",
    scaling_type="strong",
    max_nodes=8
)
orchestrator.run()
```

**What happens automatically:**
1. ✅ Code acquired (local or Git clone)
2. ✅ Dependencies detected from README
3. ✅ Code compiled with proper environment
4. ✅ Scaling tests generated and submitted
5. ✅ Efficiency reports created

➡️ **[Full Automated Workflow Guide](docs/AUTOMATED_WORKFLOW.md)**

---

### Option 2: Manual Test Definition (Advanced)

**For fine-grained control:**

```python
from pathlib import Path
from core.test_definition import Test

# Create test
test = Test(
    name="my_benchmark",
    input_file=Path("input.dat"),
    command=["./my_app", "--input", "input.dat"]
)

# Configure backend
test.set_backend(
    scheduler="slurm",
    launcher="srun",
    module_system="lmod"
)

# Configure resources
test.set_resources(
    max_nodes=64,
    procs_per_node=128,
    time_limit="02:00:00"
)

# Configure scaling
test.set_scaling(
    scaling_type="strong",
    max_nodes=64,
    initial_procs=(2, 2, 2)
)
```

## Configuration

### Test Configuration

```python
# Backend Selection
test.set_backend(
    scheduler="slurm",        # local, slurm, pbs
    launcher="srun",          # srun, mpirun, mpiexec
    module_system="lmod",     # nomod, tmod, tmod4, lmod
    build_system="cmake"      # make, cmake, autotools, easybuild, spack
)

# Resource Configuration
test.set_resources(
    max_nodes=128,
    procs_per_node=128,
    gpus_per_node=1,
    memory_per_node="100GB",
    time_limit="02:00:00",
    partition="gpu",
    account="project123"
)

# Scaling Configuration
test.set_scaling(
    scaling_type="strong",    # or "weak"
    max_nodes=128,
    initial_procs=(2, 2, 2),         # 3D decomposition
    initial_domain=(10.0, 10.0, 10.0),
    initial_cells=(256, 256, 256)
)

# Environment Setup
test.set_modules(["gcc/11.2.0", "openmpi/4.1.1"])
test.set_env({"OMP_NUM_THREADS": "1"})

# Control job submission behavior
test.set_auto_submit(True)  # Default: automatically submit jobs
# test.set_auto_submit(False)  # Only generate job scripts, don't submit
```

## Scaling Tests

### Strong Scaling

Problem size stays constant, increase parallelism:

```python
test.set_scaling(
    scaling_type="strong",
    max_nodes=64,
    initial_procs=(2, 2, 2),  # Start with 8 processes
    # Domain size stays constant
)
```

Generates configs for 1, 2, 4, 8, 16, 32, 64 nodes with:
- Increasing processor counts (alternating x/y doubling)
- Constant problem size
- Efficiency metrics calculated automatically

### Weak Scaling

Problem size grows with parallelism:

```python
test.set_scaling(
    scaling_type="weak",
    max_nodes=64,
    initial_procs=(2, 2, 2),
    # Domain and cells scale proportionally
)
```

## System Configuration

**NEW**: Define custom launchers and system configurations (similar to ReFrame):

```python
# In your leonardo_system.py
from core.registry import register_launcher, JobLauncher

@register_launcher('mpirun-mapby')
class MpirunMapbyLauncher(JobLauncher):
    def command(self, job, resource_config):
        return [
            'mpirun', '-np', str(job.num_procs),
            '--map-by', 'socket:PE=8',
            '--report-bindings'
        ]

# Define system information
site_configuration = {
    "systems": [{
        "name": "leonardo",
        "partitions": [{
            "name": "booster",
            "scheduler": "slurm",
            "launcher": "mpirun-mapby",  # Use custom launcher
            "processor": {
                "num_cpus": 32,
                "num_sockets": 1
            },
            "devices": [{
                "type": "gpu",
                "model": "A100",
                "num_devices": 4
            }]
        }]
    }],
    "environments": [{
        "name": "gcc-openmpi",
        "modules": ["gcc/11.2.0", "openmpi/4.1.1"],
        "features": ["mpi", "openmp"]
    }]
}
```

**Benefits:**
- ✅ **Auto-detection**: System detected by hostname matching
- ✅ **Validation**: Ensure procs_per_node matches actual hardware
- ✅ **Custom Launchers**: Define MPI commands with specific binding
- ✅ **Environment Management**: Organize compiler/MPI combinations

➡️ **[Full System Configuration Guide](SYSTEM_CONFIG_GUIDE.md)**

**Quick Example:**

```python
from utils.system_loader import SystemConfigLoader

# Load system configuration
loader = SystemConfigLoader(Path("leonardo_system.py"))

# Auto-create resource config from partition info
resource_config = loader.create_resource_config(
    partition_name="booster",
    max_nodes=16
)
# Automatically sets correct procs_per_node, gpus_per_node, etc.
```

## Job Submission Control

HPC-ScaleTest provides flexible control over job submission:

### Automatic Submission (Default)

Jobs are automatically submitted after script generation:

```python
test.set_auto_submit(True)  # This is the default behavior
```

### Manual Submission

Only generate job scripts, submit manually later:

```python
test.set_auto_submit(False)

# Run the test to generate job scripts
# python scaletest.py run --test my_test.py

# Later, submit the jobs manually:
# python submit_jobs.py output/my_test_strong_20251019_120000
```

### Submitting Prepared Jobs

Use the provided utility script to submit prepared jobs:

```bash
# Submit all prepared jobs in an output directory
python submit_jobs.py output/my_test_strong_20251019_120000

# Or submit a single job script
python -c "
from pathlib import Path
from utils.job_submitter import submit_single_job
job_id = submit_single_job(Path('output/my_test_strong_20251019_120000/nodes_1/job.sh'))
print(f'Submitted job: {job_id}')
"
```

## Advanced Features

### Large Job Handling

For jobs exceeding 64 nodes, the framework automatically:
- Comments out default QoS
- Enables boost QoS
- Adjusts resource allocations

### GPU Support

```python
test.set_resources(
    gpus_per_node=1,
    # Automatically adds CUDA_VISIBLE_DEVICES
)
```

### Custom Build Integration

```python
from core.factory import BackendFactory

# Build with CMake
builder = BackendFactory.create_build_system(
    BuildBackend.CMAKE,
    {'parallel_jobs': 8}
)

build_path = builder.build(
    source_dir=Path("./src"),
    flags={"CMAKE_BUILD_TYPE": "Release"}
)
```

## Backend Details

### Schedulers

- **LocalScheduler**: Runs jobs as local processes (debugging)
- **SlurmScheduler**: Full Slurm integration (sbatch, squeue, scancel, sacct)

### Launchers

- **SrunLauncher**: Slurm native launcher with GPU binding
- **MpiRunLauncher**: Generic MPI launcher (supports mpirun/mpiexec)

### Module Systems

- **NoModBackend**: No-op for systems without modules
- **TModBackend**: TCL-based modules (Environment Modules 3.x)
- **TMod4Backend**: TCL-based modules (Environment Modules 4.x)
- **LModBackend**: Lua-based Lmod with spider support

### Build Systems

- **MakeBackend**: Standard Makefile builds
- **CMakeBackend**: CMake with automatic build directory generation
- **AutotoolsBackend**: configure/make/install workflows
- **EasyBuildBackend**: EasyBuild integration with robot mode
- **SpackBackend**: Spack package manager integration

## Testing

```bash
# Run unit tests
make test

# Run specific test
python -m pytest tests/test_scaling.py -v
```

## Extending the Framework

### Adding a New Scheduler

1. Create `backends/schedulers/my_scheduler.py`
2. Inherit from `AbstractScheduler`
3. Implement required methods: `submit`, `monitor`, `cancel`, `wait`
4. Add to `SchedulerBackend` enum in `core/types.py`
5. Register in `BackendFactory.create_scheduler()`

Example:

```python
from core.abstracts import AbstractScheduler
from core.types import JobID, JobStatus

class MyScheduler(AbstractScheduler):
    def submit(self, job_script: Path) -> JobID:
        # Submit job using your scheduler's CLI
        pass
    
    def monitor(self, job_id: JobID) -> JobStatus:
        # Check job status
        pass
    
    def cancel(self, job_id: JobID) -> bool:
        # Cancel job
        pass
    
    def wait(self, job_id: JobID, timeout: int = None) -> JobResults:
        # Wait for completion and gather results
        pass
```

### Adding a New Build System

Similar process for launchers, module systems, and build systems.