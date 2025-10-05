# HPC-ScaleTest Usage Guide

Complete guide with practical examples for using HPC-ScaleTest.

## Table of Contents

1. [Basic Workflow](#basic-workflow)
2. [Test Definition Patterns](#test-definition-patterns)
3. [Backend Configuration](#backend-configuration)
4. [Scaling Studies](#scaling-studies)
5. [Result Analysis](#result-analysis)
6. [Advanced Scenarios](#advanced-scenarios)

## Basic Workflow

### Step 1: Create Test Definition

```python
# my_first_test.py
from pathlib import Path
from core.test_definition import Test

# Create test instance
test = Test(
    name="my_benchmark",
    input_file=Path("./inputs/config.dat"),
    command=["./my_application", "--mode", "fast"]
)

# Configure for local testing
test.set_backend(scheduler="local", launcher="mpirun")
test.set_resources(max_nodes=2, procs_per_node=4)
test.set_scaling(
    scaling_type="strong",
    max_nodes=2,
    initial_procs=(1, 1, 1)
)
```

### Step 2: Run Test

```bash
python scaletest.py run --test my_first_test.py
```

### Step 3: Examine Results

```bash
cd output/my_benchmark_strong_*/
cat summary.json
cat efficiency_report.txt
```

## Test Definition Patterns

### Pattern 1: Inline Configuration

```python
from pathlib import Path
from core.test_definition import Test

# All configuration in one place
test = Test("inline_test", Path("input.dat"), ["./app"]) \
    .set_backend(scheduler="slurm", launcher="srun") \
    .set_resources(max_nodes=32, procs_per_node=128) \
    .set_scaling(scaling_type="strong", max_nodes=32,
                 initial_procs=(2, 2, 2))
```

### Pattern 2: Function Factory

```python
def create_cfd_test(max_nodes=64):
    """Factory function for CFD tests."""
    test = Test("cfd_test", Path("cfd.inp"), ["./cfd_solver"])
    
    test.set_modules(["gcc/11.2", "openmpi/4.1", "hdf5/1.12"])
    test.set_env({
        "OMP_NUM_THREADS": "1",
        "HDF5_USE_FILE_LOCKING": "FALSE"
    })
    
    test.set_backend(
        scheduler="slurm",
        launcher="srun",
        module_system="lmod"
    )
    
    test.set_resources(
        max_nodes=max_nodes,
        procs_per_node=128,
        time_limit="04:00:00",
        partition="compute",
        account="CFD_PROJECT"
    )
    
    test.set_scaling(
        scaling_type="strong",
        max_nodes=max_nodes,
        initial_procs=(4, 4, 4),
        initial_domain=(100.0, 100.0, 100.0),
        initial_cells=(512, 512, 512)
    )
    
    return test

# Create test
test = create_cfd_test(max_nodes=128)
```

### Pattern 3: Class-Based

```python
class MyBenchmarkTest(Test):
    """Custom test class with preset configuration."""
    
    def __init__(self, max_nodes=64):
        super().__init__(
            name="custom_benchmark",
            input_file=Path("input.cfg"),
            command=["./benchmark_app"]
        )
        
        self._configure(max_nodes)
    
    def _configure(self, max_nodes):
        """Apply standard configuration."""
        self.set_modules(["intel/2021", "impi/2021"])
        self.set_backend(scheduler="slurm", launcher="srun")
        self.set_resources(
            max_nodes=max_nodes,
            procs_per_node=128,
            time_limit="02:00:00"
        )
        self.set_scaling(
            scaling_type="weak",
            max_nodes=max_nodes,
            initial_procs=(2, 2, 2),
            initial_domain=(10, 10, 10),
            initial_cells=(256, 256, 256)
        )

# Usage
test = MyBenchmarkTest(max_nodes=256)
```

## Backend Configuration

### Local Development

```python
# For debugging on workstation
test.set_backend(
    scheduler="local",      # Run as local processes
    launcher="mpirun",      # Use mpirun
    module_system="nomod",  # No modules
    build_system="make"     # Simple make
)

test.set_resources(
    max_nodes=1,
    procs_per_node=4  # Match your core count
)
```

### Slurm Cluster

```python
# Production HPC system
test.set_backend(
    scheduler="slurm",
    launcher="srun",
    module_system="lmod",
    build_system="cmake"
)

test.set_resources(
    max_nodes=256,
    procs_per_node=128,
    time_limit="08:00:00",
    partition="standard",
    account="project_123",
    memory_per_node="256GB"
)
```

### GPU Cluster

```python
# GPU-accelerated application
test.set_backend(
    scheduler="slurm",
    launcher="srun",
    module_system="lmod"
)

test.set_resources(
    max_nodes=16,
    procs_per_node=4,
    gpus_per_node=4,  # 4 GPUs per node
    time_limit="12:00:00",
    partition="gpu",
    account="gpu_project"
)

test.set_modules([
    "gcc/11.2",
    "cuda/11.7",
    "openmpi/4.1+cuda"
])

test.set_env({
    "CUDA_VISIBLE_DEVICES": "0,1,2,3",
    "CUDA_DEVICE_ORDER": "PCI_BUS_ID"
})
```

## Scaling Studies

### Strong Scaling Example

Problem: Fixed-size 3D CFD simulation

```python
test = Test("cfd_strong", Path("cfd_512cubed.inp"), ["./cfd"])

test.set_scaling(
    scaling_type="strong",
    max_nodes=128,
    # Start with 8 processes (2x2x2)
    initial_procs=(2, 2, 2),
    # Fixed problem size: 512^3 cells in 100x100x100 domain
    initial_domain=(100.0, 100.0, 100.0),
    initial_cells=(512, 512, 512)
)

# Generates: 1, 2, 4, 8, 16, 32, 64, 128 nodes
# Process decomposition grows: 2x2x2 -> 4x2x2 -> 4x4x2 -> 8x4x2 -> ...
# Problem size stays constant
```

Expected behavior:
- Runtime should decrease as nodes increase
- Efficiency = (T1 / (N * TN)) where T1 is baseline time, N is node count
- Perfect scaling: 2x nodes = 0.5x time

### Weak Scaling Example

Problem: Maintaining constant work per processor

```python
test = Test("molecular_weak", Path("md.inp"), ["./md_sim"])

test.set_scaling(
    scaling_type="weak",
    max_nodes=256,
    # Start with 8 processes
    initial_procs=(2, 2, 2),
    # Problem grows with processor count
    initial_domain=(50.0, 50.0, 50.0),
    initial_cells=(256, 256, 256)
)

# Generates configurations where:
# - 1 node: 2x2x2 procs, 50x50x50 domain, 256^3 cells
# - 2 nodes: 4x2x2 procs, 100x50x50 domain, 512x256x256 cells
# - 4 nodes: 4x4x2 procs, 100x100x50 domain, 512x512x256 cells
# Problem size grows proportionally
```

Expected behavior:
- Runtime should stay approximately constant
- Efficiency = T1 / TN (ratio of baseline to current time)
- Perfect scaling: constant time regardless of node count

### Mixed CPU/GPU Scaling

```python
test = Test("hybrid_app", Path("input.dat"), ["./hybrid"])

test.set_resources(
    max_nodes=64,
    procs_per_node=32,  # CPU processes
    gpus_per_node=2,    # GPU accelerators
)

test.set_scaling(
    scaling_type="strong",
    max_nodes=64,
    initial_procs=(2, 2, 1),  # Adjust for GPU topology
    initial_domain=(200, 200, 100),
    initial_cells=(1024, 1024, 512)
)

test.set_env({
    "OMP_NUM_THREADS": "8",  # CPU threads per MPI rank
    "CUDA_VISIBLE_DEVICES": "0,1"
})
```

## Result Analysis

### Examining Output Structure

```bash
output/my_test_strong_20250104_143022/
├── nodes_1/
│   ├── job.sh              # Generated job script
│   ├── slurm-12345.out     # Job stdout
│   ├── slurm-12345.err     # Job stderr
│   ├── results.json        # Parsed results
│   └── input.dat           # Copied input file
├── nodes_2/
├── nodes_4/
├── ...
├── summary.json            # Aggregated results
└── efficiency_report.txt   # Performance analysis
```

### Parsing Results Programmatically

```python
import json
from pathlib import Path

# Load summary
result_dir = Path("output/my_test_strong_20250104_143022")
with open(result_dir / "summary.json") as f:
    summary = json.load(f)

# Analyze results
for n_nodes, result in summary['results'].items():
    if result['status'] == 'COMPLETED':
        runtime = result['runtime']
        procs = result['total_procs']
        print(f"{n_nodes} nodes ({procs} procs): {runtime:.2f}s")

# Calculate custom metrics
baseline = summary['results']['1']
baseline_time = baseline['runtime']

for n_nodes, result in summary['results'].items():
    if result['status'] == 'COMPLETED':
        n = int(n_nodes)
        speedup = baseline_time / result['runtime']
        efficiency = (speedup / n) * 100
        print(f"Node {n_nodes}: Speedup={speedup:.2f}x, Efficiency={efficiency:.1f}%")
```

### Plotting Results

```python
import json
import matplotlib.pyplot as plt
from pathlib import Path

def plot_scaling(result_dir):
    """Plot scaling results."""
    with open(result_dir / "summary.json") as f:
        summary = json.load(f)
    
    nodes = []
    times = []
    
    for n_nodes, result in sorted(summary['results'].items(), key=lambda x: int(x[0])):
        if result['status'] == 'COMPLETED':
            nodes.append(int(n_nodes))
            times.append(result['runtime'])
    
    # Plot runtime vs nodes
    plt.figure(figsize=(10, 6))
    plt.plot(nodes, times, 'o-', linewidth=2, markersize=8)
    plt.xlabel('Number of Nodes')
    plt.ylabel('Runtime (seconds)')
    plt.title('Strong Scaling Performance')
    plt.grid(True)
    plt.xscale('log', base=2)
    plt.savefig(result_dir / 'scaling_plot.png', dpi=300)
    
    # Plot speedup
    baseline_time = times[0]
    speedups = [baseline_time / t for t in times]
    ideal_speedup = nodes
    
    plt.figure(figsize=(10, 6))
    plt.plot(nodes, speedups, 'o-', label='Actual', linewidth=2, markersize=8)
    plt.plot(nodes, ideal_speedup, '--', label='Ideal', linewidth=2)
    plt.xlabel('Number of Nodes')
    plt.ylabel('Speedup')
    plt.title('Speedup Curve')
    plt.legend()
    plt.grid(True)
    plt.xscale('log', base=2)
    plt.yscale('log', base=2)
    plt.savefig(result_dir / 'speedup_plot.png', dpi=300)

# Usage
plot_scaling(Path("output/my_test_strong_20250104_143022"))
```

## Advanced Scenarios

### Scenario 1: Parametric Study

Run same test with different input parameters:

```python
# parametric_study.py
from pathlib import Path
from core.test_definition import Test

def create_parametric_test(param_value):
    """Create test for specific parameter value."""
    test = Test(
        name=f"param_study_{param_value}",
        input_file=Path(f"inputs/config_{param_value}.dat"),
        command=["./app", "--param", str(param_value)]
    )
    
    test.set_backend(scheduler="slurm", launcher="srun")
    test.set_resources(max_nodes=32, procs_per_node=128)
    test.set_scaling(scaling_type="strong", max_nodes=32,
                     initial_procs=(2, 2, 2))
    
    return test

# Create tests for different parameters
for param in [0.1, 0.5, 1.0, 2.0, 5.0]:
    test = create_parametric_test(param)
    # Save or run each test
```

### Scenario 2: Multi-Stage Pipeline

Build, then test:

```python
from pathlib import Path
from core.test_definition import Test
from core.factory import BackendFactory
from core.types import BuildBackend

# Stage 1: Build application
builder = BackendFactory.create_build_system(
    BuildBackend.CMAKE,
    options={'parallel_jobs': 16}
)

build_path = builder.build(
    source_dir=Path("./src"),
    flags={
        "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_CXX_COMPILER": "mpicxx",
        "ENABLE_GPU": "ON"
    }
)

install_path = builder.install(build_path, Path("./install"))

# Stage 2: Run tests with built executable
test = Test(
    name="pipeline_test",
    input_file=Path("input.dat"),
    command=[str(install_path / "bin" / "app")]
)

test.set_resources(max_nodes=64, procs_per_node=128)
test.set_scaling(scaling_type="strong", max_nodes=64,
                 initial_procs=(4, 4, 4))
```

### Scenario 3: Conditional Configuration

Adapt based on environment:

```python
import os
from pathlib import Path
from core.test_definition import Test

def create_adaptive_test():
    """Create test that adapts to environment."""
    test = Test("adaptive", Path("input.dat"), ["./app"])
    
    # Detect if running on HPC or local
    if os.environ.get('SLURM_JOB_ID'):
        # Running on Slurm cluster
        test.set_backend(scheduler="slurm", launcher="srun",
                        module_system="lmod")
        test.set_resources(max_nodes=128, procs_per_node=128,
                          time_limit="04:00:00", partition="compute")
        test.set_modules(["gcc/11.2", "openmpi/4.1"])
    else:
        # Running locally
        test.set_backend(scheduler="local", launcher="mpirun",
                        module_system="nomod")
        test.set_resources(max_nodes=1, procs_per_node=4)
    
    test.set_scaling(scaling_type="strong", max_nodes=test.resource_config.max_nodes,
                     initial_procs=(2, 2, 2))
    
    return test

test = create_adaptive_test()
```

### Scenario 4: Custom Job Script Templates

Override job script generation:

```python
from pathlib import Path
from core.test_definition import Test

test = Test("custom_script", Path("input.dat"), ["./app"])

# Standard configuration
test.set_backend(scheduler="slurm", launcher="srun")
test.set_resources(max_nodes=32, procs_per_node=128)

# Custom environment setup
test.set_env({
    # Custom module paths
    "MODULEPATH": "/custom/modules:$MODULEPATH",
    
    # Performance tuning
    "OMP_NUM_THREADS": "1",
    "MPI_BUFFER_SIZE": "32768",
    "UCX_TLS": "rc_x,self",
    
    # I/O optimization
    "ROMIO_HINTS": "./romio_hints.txt",
    "HDF5_USE_FILE_LOCKING": "FALSE",
    
    # Debug options
    "MPICH_DBG": "class,file",
})

# Custom modules with versions
test.set_modules([
    "gcc/11.2.0",
    "openmpi/4.1.1+ucx",
    "hdf5/1.12.1+parallel",
    "fftw/3.3.10",
    "boost/1.77.0"
])
```

### Scenario 5: Convergence Study

Test with increasing resolution:

```python
def create_convergence_test(resolution):
    """Create test with specific grid resolution."""
    cells = 128 * resolution  # 128, 256, 512, 1024...
    
    test = Test(
        name=f"convergence_res{resolution}",
        input_file=Path(f"inputs/grid_{cells}.dat"),
        command=["./solver", "--grid-size", str(cells)]
    )
    
    # Scale resources with problem size
    base_nodes = 4
    nodes = base_nodes * (resolution ** 2)  # More nodes for finer grids
    
    test.set_backend(scheduler="slurm", launcher="srun")
    test.set_resources(
        max_nodes=min(nodes, 256),  # Cap at 256 nodes
        procs_per_node=128,
        time_limit="08:00:00"
    )
    
    test.set_scaling(
        scaling_type="strong",
        max_nodes=min(nodes, 256),
        initial_procs=(4, 4, 4),
        initial_domain=(100, 100, 100),
        initial_cells=(cells, cells, cells)
    )
    
    return test

# Run convergence study
for res in [1, 2, 4, 8]:
    test = create_convergence_test(res)
    # Run each test
```

### Scenario 6: Checkpointing and Restart

Handle long-running jobs with checkpoints:

```python
test = Test("checkpoint_test", Path("input.dat"), ["./app"])

# Set checkpoint parameters
test.set_env({
    "CHECKPOINT_INTERVAL": "3600",  # Every hour
    "CHECKPOINT_DIR": "./checkpoints",
    "RESTART_FROM_CHECKPOINT": "false"
})

# Extended time with checkpointing
test.set_resources(
    max_nodes=64,
    procs_per_node=128,
    time_limit="24:00:00"  # Long job
)

# For restart runs, modify:
# test.set_env({
#     "RESTART_FROM_CHECKPOINT": "true",
#     "CHECKPOINT_FILE": "./checkpoints/latest.chk"
# })
```

## Best Practices

### 1. Version Control Test Definitions

```python
# test_v1.0.py
"""
Test configuration for MyApp v1.0
Date: 2025-01-04
Author: Your Name
Description: Baseline strong scaling test
"""

from pathlib import Path
from core.test_definition import Test

VERSION = "1.0.0"
APP_VERSION = "2.5.3"

test = Test(f"myapp_v{VERSION}", Path("input.dat"), ["./myapp"])
# ... configuration ...
```

### 2. Document Expected Performance

```python
def create_baseline_test():
    """
    Baseline performance test.
    
    Expected results (as of 2025-01-04):
    - 1 node: ~100s
    - 16 nodes: ~7s (93% efficiency)
    - 64 nodes: ~2.5s (62% efficiency)
    
    Known issues:
    - Communication overhead increases above 64 nodes
    - May need MPI tuning for > 128 nodes
    """
    test = Test("baseline", Path("input.dat"), ["./app"])
    # ... configuration ...
    return test
```

### 3. Validate Before Running

Always validate before large runs:

```bash
# Validate test definition
python scaletest.py validate --test my_test.py

# Test locally first
python scaletest.py run --test my_test.py --backend local --max-nodes 1

# Then scale up
python scaletest.py run --test my_test.py --max-nodes 128
```

### 4. Monitor Resource Usage

Add monitoring to job scripts:

```python
test.set_env({
    # Enable Slurm profiling
    "SLURM_PROFILE": "All",
    
    # Monitor memory
    "MALLOC_TRIM_THRESHOLD_": "131072",
    
    # Track MPI stats
    "MPICH_ENV_DISPLAY": "1"
})
```

### 5. Archive Results

```bash
# After test run
cd output/
tar -czf my_test_results_$(date +%Y%m%d).tar.gz my_test_strong_*/
mv my_test_results_*.tar.gz ~/archive/
```

## Troubleshooting Tips

### Common Issues and Solutions

**Issue**: Jobs fail immediately  
**Solution**: Check job script permissions and paths
```bash
ls -la output/test_*/nodes_*/job.sh
cat output/test_*/nodes_1/slurm-*.err
```

**Issue**: Module load failures  
**Solution**: Verify modules exist
```bash
module avail gcc
module spider openmpi
```

**Issue**: Inconsistent timing  
**Solution**: Add warm-up runs or exclude first run
```python
# Run benchmark twice, use second timing
command=["./app", "--warmup", "1", "--runs", "2"]
```

**Issue**: Poor scaling  
**Solution**: Check for:
- I/O bottlenecks (use profilers)
- Load imbalance (check per-rank times)
- Communication overhead (MPI profiling)

## Summary

This guide covered:
- ✅ Creating test definitions
- ✅ Configuring backends for different systems
- ✅ Running strong and weak scaling studies
- ✅ Analyzing and plotting results
- ✅ Advanced scenarios and best practices

For more information, see README.md and API documentation.
