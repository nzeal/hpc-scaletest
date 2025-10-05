# HPC-ScaleTest

A modular Python framework for running benchmark scaling tests on heterogeneous HPC systems with support for CPU and GPU resources.

## Features

- **Layered Architecture**: Clean separation between user-facing API and backend implementations
- **Scheduler Agnostic**: Write tests once, run on any system (Local, Slurm, PBS, etc.)
- **Pluggable Backends**: Easily swap schedulers, launchers, module systems, and build tools
- **Scaling Tests**: Built-in support for strong and weak scaling studies
- **Heterogeneous Support**: Handle mixed CPU/GPU resources
- **Automatic Job Generation**: Generate and submit jobs with proper resource allocations
- **Result Aggregation**: Collect and analyze performance metrics

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

### 1. Define a Test

Create `my_test.py`:

```python
from pathlib import Path
from core.test_definition import Test

def create_test():
    test = Test(
        name="my_benchmark",
        input_file=Path("./input.dat"),
        command=["./my_app"]
    )
    
    test.set_backend(
        scheduler="slurm",
        launcher="srun"
    )
    
    test.set_resources(
        max_nodes=64,
        procs_per_node=128,
        time_limit="01:00:00"
    )
    
    test.set_scaling(
        scaling_type="strong",
        max_nodes=64,
        initial_procs=(2, 2, 2),
        initial_domain=(10.0, 10.0, 10.0),
        initial_cells=(256, 256, 256)
    )
    
    return test

test = create_test()
```

### 2. Run the Test

```bash
# Run with Slurm
python scaletest.py run --test my_test.py --scaling strong --max-nodes 64

# Run locally for debugging
python scaletest.py run --test my_test.py --backend local --max-nodes 4

# Validate test definition
python scaletest.py validate --test my_test.py

# List available backends
python scaletest.py list-backends
```

### 3. Check Results

Results are organized in timestamped directories:

```
output/
└── my_benchmark_strong_20250104_143022/
    ├── nodes_1/
    │   ├── job.sh
    │   ├── slurm-12345.out
    │   └── results.json
    ├── nodes_2/
    ├── nodes_4/
    ├── ...
    ├── summary.json
    └── efficiency_report.txt
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

Run unit tests:

```bash
# Test scaling engine
python -m unittest tests.test_scaling

# Test backends
python -m unittest tests.test_backends

# Run all tests
python -m unittest discover tests
```

## Project Structure

```
hpc_scaletest/
├── scaletest.py              # Main CLI entry point
├── core/
│   ├── abstracts.py          # Abstract base classes
│   ├── types.py              # Type definitions and enums
│   ├── config.py             # Configuration dataclasses
│   ├── factory.py            # Backend factory
│   └── test_definition.py    # User-facing Test class
├── backends/
│   ├── schedulers/
│   │   ├── local.py
│   │   └── slurm.py
│   ├── launchers/
│   │   ├── srun.py
│   │   └── mpirun.py
│   ├── modules/
│   │   ├── nomod.py
│   │   ├── tmod.py
│   │   ├── tmod4.py
│   │   └── lmod.py
│   └── builds/
│       ├── make.py
│       ├── cmake.py
│       ├── autotools.py
│       ├── easybuild.py
│       └── spack.py
├── engine/
│   ├── scaling.py            # Scaling configuration generator
│   └── runner.py             # Test execution engine
├── utils/
│   ├── file_utils.py         # File manipulation utilities
│   └── logging_config.py     # Logging setup
├── tests/
│   ├── test_scaling.py
│   └── test_backends.py
└── examples/
    ├── example_test.py
    └── config.yaml
```

## CLI Reference

### Commands

```bash
# Run a test
scaletest.py run --test <file> [options]

# Validate test definition
scaletest.py validate --test <file>

# List available backends
scaletest.py list-backends
```

### Options

```
--test <file>           Test definition file (required)
--scaling <type>        Scaling type: strong or weak
--max-nodes <n>         Maximum nodes to scale to
--backend <name>        Scheduler backend
--output <dir>          Output directory
-v, --verbose           Enable verbose logging
```

## Examples

### Example 1: Local Testing

```python
from pathlib import Path
from core.test_definition import Test

test = Test("simple_test", Path("input.dat"), ["./app"])
test.set_backend(scheduler="local", launcher="mpirun")
test.set_resources(max_nodes=2, procs_per_node=4)
test.set_scaling(scaling_type="strong", max_nodes=2, 
                 initial_procs=(1, 1, 1))
```

Run: `python scaletest.py run --test simple_test.py`

### Example 2: Slurm Cluster

```python
test = Test("hpc_test", Path("input.dat"), ["./hpc_app"])
test.set_modules(["gcc/11.2", "openmpi/4.1"])
test.set_backend(scheduler="slurm", launcher="srun", 
                 module_system="lmod")
test.set_resources(max_nodes=128, procs_per_node=128, 
                   time_limit="04:00:00", partition="compute")
test.set_scaling(scaling_type="weak", max_nodes=128,
                 initial_procs=(2, 2, 2),
                 initial_domain=(10, 10, 10),
                 initial_cells=(256, 256, 256))
```

Run: `python scaletest.py run --test hpc_test.py --max-nodes 128`

### Example 3: GPU Application

```python
test = Test("gpu_benchmark", Path("input.dat"), ["./gpu_app"])
test.set_modules(["cuda/11.7", "openmpi/4.1"])
test.set_backend(scheduler="slurm", launcher="srun")
test.set_resources(max_nodes=16, procs_per_node=4, 
                   gpus_per_node=4, partition="gpu")
test.set_scaling(scaling_type="strong", max_nodes=16,
                 initial_procs=(1, 1, 1))
test.set_env({"CUDA_VISIBLE_DEVICES": "0,1,2,3"})
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

## Performance Analysis

The framework generates efficiency reports for strong scaling tests:

```
Strong Scaling Efficiency Report
==================================================

Nodes      Procs      Time(s)      Speedup      Efficiency
------------------------------------------------------------
1          128        100.00       1.00         100.0%
2          256        52.00        1.92         96.2%
4          512        28.00        3.57         89.3%
8          1024       16.00        6.25         78.1%
16         2048       10.00        10.00        62.5%
```

Access results programmatically:

```python
import json
from pathlib import Path

# Load summary
with open('output/test_strong_*/summary.json') as f:
    summary = json.load(f)

# Process results
for n_nodes, result in summary['results'].items():
    print(f"{n_nodes} nodes: {result['runtime']}s")
```

## Troubleshooting

### Common Issues

1. **Module not found**: Ensure modules are available on your system
   ```bash
   module avail  # Check available modules
   ```

2. **Job submission fails**: Check scheduler configuration
   ```bash
   scaletest.py validate --test mytest.py
   ```

3. **Permission denied**: Make sure scripts are executable
   ```bash
   chmod +x scaletest.py
   ```

4. **Import errors**: Ensure project structure is correct and Python path is set

### Debug Mode

Enable verbose logging:

```bash
python scaletest.py run --test mytest.py --verbose
```

## Contributing

1. Follow PEP 8 style guidelines
2. Add docstrings to all classes and methods
3. Write unit tests for new features
4. Update documentation

## License

[Your License Here]

## Contact

[Your Contact Information]

## Acknowledgments

This framework was designed for HPC scaling studies with emphasis on:
- Portability across different HPC systems
- Ease of use for researchers
- Extensibility for new backends
- Reproducible performance testing
