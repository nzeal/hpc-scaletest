# YAML Configuration Guide for HPC-ScaleTest

## Overview

HPC-ScaleTest now supports YAML configuration files, providing a cleaner and more maintainable alternative to long command-line arguments. This guide explains how to use YAML configuration for automated scaling tests.

## Quick Start

### 1. Create a Configuration File

Create a file named `run.yaml` (or any name you prefer):

```yaml
repository: https://github.com/thenitinshukla/iPIC3D-CPU-NS.git
scaling: weak
nodes: 4
partition: dcgp
account: cin_staff
modules:
  - hdf5/1.14.3--intel-oneapi-mpi--2021.12.1--oneapi--2024.1.0
input_file: inputfiles/os-stdin
verbose: true
```

### 2. Run the Test

```bash
python hpc_auto.py --config run.yaml
```

That's it! The framework will:
- Clone the repository (if not already present)
- Auto-detect system resources (CPUs, GPUs, memory)
- Detect build system and dependencies
- Compile the code
- Generate scaling configurations
- Submit jobs and collect results

## Key Features

### 🎯 Dynamic System Detection

The framework automatically detects:
- **Total nodes** available in your HPC system
- **Cores per node** (from Slurm, PBS, or `/proc/cpuinfo`)
- **GPUs per node** (via `nvidia-smi` or scheduler info)
- **Memory per node** (from scheduler or `/proc/meminfo`)
- **Scheduler type** (Slurm, PBS, or local)

You only need to specify what you want to override!

### 📊 Intelligent Scaling Configuration

#### Three Operating Modes

HPC-ScaleTest supports three distinct operating modes:

##### 1. **Baseline Mode** (No Scaling - Single Job)
When no scaling factors are defined, run a single simulation with exact initial parameters:
```yaml
repository: https://github.com/user/repo.git
partition: dcgp
account: cin_staff

# Initial parameters - used EXACTLY as provided
initial_procs: [14, 8, 1]  # 112 MPI processes
initial_domain: [84.0, 48.0, 1.0]
initial_cells: [840, 480, 1]
particles_per_cell: {x: 20, y: 20, z: 20}

# No scaling factors = single job with exact parameters
```

**Output**: Single job folder (e.g., `output/repo_initial/`) with unmodified parameters.

##### 2. **Weak Scaling Mode** (Problem Size Grows with Resources)
Problem size scales with number of resources to maintain constant work per process:

**Using Scaling Factors** (Recommended):
```yaml
scaling: weak
weak_scaling_factors: [1, 2, 4, 8]  # Explicit scaling factors

# Initial parameters (used exactly for S=1)
initial_procs: [14, 8, 1]
initial_domain: [84.0, 48.0, 1.0]
initial_cells: [840, 480, 1]
particles_per_cell: {x: 20, y: 20, z: 20}  # CONSTANT across all factors
```

**Rules**:
- **S=1**: Uses EXACT initial values (no modification)
- **S>1**: Scales domain, grid, and MPI decomposition by ∛S (cubic root)
- **Particles per cell**: Remains CONSTANT across all scaling factors
- **Formula**: `scaled_value = initial_value × ∛S`

**Example for S=2**:
- Decomposition: [14×∛2, 8×∛2, 1×∛2] ≈ [18, 10, 1] = 224 procs
- Domain: [84.0×∛2, 48.0×∛2, 1.0×∛2] ≈ [105.8, 60.5, 1.26]
- Cells: [840×∛2, 480×∛2, 1×∛2] ≈ [1058, 605, 1]
- Particles/cell: [20, 20, 20] ✅ UNCHANGED

**Using Node Sequence** (Legacy):
```yaml
scaling: weak
nodes: 4  # Generates node sequence: 1, 2, 4

initial_procs: [2, 2, 2]
initial_domain: [10.0, 10.0, 10.0]
initial_cells: [256, 256, 256]
```

**Output**: Folders `weak_1/`, `weak_2/`, `weak_4/`, `weak_8/` (or `weak_nodes_1/`, etc. for node-based)

##### 3. **Strong Scaling Mode** (Fixed Problem Size)
Problem size remains CONSTANT while increasing computational resources:

**Using Scaling Factors** (Recommended):
```yaml
scaling: strong
strong_scaling_factors: [1, 2, 4, 8]

# Initial parameters
initial_procs: [14, 8, 1]
initial_domain: [84.0, 48.0, 1.0]  # CONSTANT
initial_cells: [840, 480, 1]        # CONSTANT
particles_per_cell: {x: 20, y: 20, z: 20}  # CONSTANT
```

**Rules**:
- **Domain**: Remains CONSTANT for all scaling factors
- **Grid**: Remains CONSTANT for all scaling factors
- **Particles per cell**: Remains CONSTANT
- **MPI decomposition**: Scales by ∛S
- **Nodes**: Auto-computed from total processes

**Example for S=4**:
- Decomposition: [14×∛4, 8×∛4, 1×∛4] ≈ [22, 13, 2] = 448 procs
- Domain: [84.0, 48.0, 1.0] ✅ UNCHANGED
- Cells: [840, 480, 1] ✅ UNCHANGED
- Particles/cell: [20, 20, 20] ✅ UNCHANGED

**Using Node Sequence** (Legacy):
```yaml
scaling: strong
nodes: 16
initial_procs: [2, 2, 2]
initial_cells: [512, 512, 512]
```

**Output**: Folders `strong_1/`, `strong_2/`, `strong_4/`, `strong_8/`

---

#### Scaling Configuration Parameters

```yaml
# Method 1: Using explicit scaling factors (RECOMMENDED)
weak_scaling_factors: [1, 2, 4, 8]      # For weak scaling
strong_scaling_factors: [1, 2, 4, 8]    # For strong scaling

# Method 2: Using node sequence (legacy)
scaling: weak | strong
nodes: 16  # Generates power-of-2 sequence: 1, 2, 4, 8, 16

# Initial parameters (REQUIRED)
initial_procs: [px, py, pz]             # MPI decomposition
initial_domain: [Lx, Ly, Lz]            # Physical domain size
initial_cells: [nx, ny, nz]             # Grid resolution
particles_per_cell: {x: npx, y: npy, z: npz}  # Particles per cell (CONSTANT)
```

#### Key Principles

1. **Initial values are sacred**: For S=1 (or first node), initial parameters are used EXACTLY as provided
2. **Particles per cell are constant**: Never modified in any scaling mode
3. **Cubic root scaling for 3D**: Maintains proper aspect ratios and work distribution
4. **Validation**: All scaling factors must be ≥ 1; initial parameters must be positive

The framework applies **cubic root scaling** for 3D problems to maintain proper work-per-process ratios.

### 🔧 Automatic Code Analysis

The framework analyzes your repository to detect:
- Build system (CMake, Make, Autotools)
- Required modules (MPI, HDF5, CUDA, etc.)
- Compilers needed
- Dependencies

## Configuration Reference

### Required Fields

```yaml
# Source code location (required)
repository: https://github.com/user/repo.git
# OR
repository: /path/to/local/code
```

### Scaling Configuration

```yaml
# Scaling type (default: strong)
scaling: weak | strong

# Maximum nodes (default: 4)
nodes: 16

# Initial processor decomposition
initial_procs: [2, 2, 2]

# For weak scaling - initial domain size
initial_domain: [10.0, 10.0, 10.0]

# For weak scaling - initial grid resolution
initial_cells: [256, 256, 256]
```

### Hardware Configuration

```yaml
# Hardware type (default: cpu)
hardware: cpu | gpu

# Override auto-detected values
procs_per_node: 128
gpus_per_node: 4
```

### Resource Configuration

```yaml
# Required for job submission
partition: dcgp
account: cin_staff

# Optional
time_limit: "02:00:00"
```

### Environment & Modules

```yaml
# Modules to load
modules:
  - gcc/11.2.0
  - openmpi/4.1.1
  - hdf5/1.14.3
  - cuda/11.8

# Module system (auto-detected)
module_system: lmod | tmod | tmod4 | nomod
```

### Build Configuration

```yaml
# Build system (auto-detected if not specified)
build_system: cmake | make | autotools

# Skip building (use existing executable)
no_build: false
```

### Input Files

```yaml
# Path to input file or directory
input_file: inputfiles/os-stdin

# Specific input file name (if directory)
input_name: os-stdin

# Disable auto-detection
no_auto_input: false
```

### Scheduler Configuration

```yaml
# Scheduler type (auto-detected)
scheduler: slurm | pbs | local

# MPI launcher (auto-selected)
launcher: srun | mpirun | mpiexec
```

### Behavior Configuration

```yaml
# Logging
verbose: true

# Job submission
auto_submit: true

# Report generation
no_reports: false

# Cleanup after completion
cleanup: false

# Custom test name
test_name: my_scaling_test
```

### Output Configuration

```yaml
# Output directories
output_dir: output
workspace_dir: workspace
```

## Advanced Usage

### Overriding YAML with CLI

Command-line arguments override YAML configuration:

```bash
python hpc_auto.py --config run.yaml --nodes 32 --partition gpu_prod
```

### Using System Configuration Files

For site-specific settings:

```yaml
system_config: leonardo_system.py
partition_name: booster
use_custom_launcher: mpirun-mapby
```

### Mixed Configuration

You can use both YAML and CLI for different purposes:

```bash
# Use YAML for standard settings, CLI for one-off changes
python hpc_auto.py --config run.yaml --scaling strong --nodes 8
```

## Complete Example

Here's a full configuration for a weak scaling test:

```yaml
# run.yaml - Weak Scaling Test for iPIC3D

# Source
repository: https://github.com/thenitinshukla/iPIC3D-CPU-NS.git

# Scaling Configuration
scaling: weak
nodes: 16
initial_procs: [2, 2, 2]
initial_domain: [10.0, 10.0, 10.0]
initial_cells: [256, 256, 256]

# Resources
partition: dcgp
account: cin_staff
time_limit: "04:00:00"

# Environment
modules:
  - hdf5/1.14.3--intel-oneapi-mpi--2021.12.1--oneapi--2024.1.0

# Input
input_file: inputfiles/os-stdin

# Behavior
verbose: true
auto_submit: true
```

Run with:
```bash
python hpc_auto.py --config run.yaml
```

## Comparison: Old vs New

### Old Way (Long CLI)
```bash
python hpc_auto.py https://github.com/user/repo.git \
  --scaling weak --nodes 16 \
  --initial-procs 2,2,2 \
  --initial-domain 10.0,10.0,10.0 \
  --initial-cells 256,256,256 \
  --partition dcgp --account cin_staff \
  --modules 'hdf5/1.14.3--intel-oneapi-mpi--2021.12.1--oneapi--2024.1.0' \
  --input-file inputfiles/os-stdin \
  --verbose
```

### New Way (YAML)
```bash
python hpc_auto.py --config run.yaml
```

Where `run.yaml` contains all the parameters in a readable format!

## Tips and Best Practices

1. **Start with minimal config**: Only specify repository, partition, and account. Let auto-detection handle the rest.

2. **Use comments**: YAML supports comments - document your choices!

3. **Version control**: Keep your `run.yaml` files in Git to track test configurations.

4. **Template approach**: Create template YAML files for different scenarios:
   - `run_strong.yaml` - Strong scaling template
   - `run_weak.yaml` - Weak scaling template  
   - `run_gpu.yaml` - GPU benchmark template

5. **Combine with CLI**: Use YAML for standard settings, CLI for quick overrides.

## Troubleshooting

### YAML Parse Errors
```bash
# Check YAML syntax
python -c "import yaml; yaml.safe_load(open('run.yaml'))"
```

### Missing Repository
```yaml
# Error: No source specified
# Solution: Add repository field
repository: /path/to/code
```

### Auto-detection Not Working
```yaml
# Override auto-detected values
procs_per_node: 128
scheduler: slurm
```

## Migration from CLI

To convert existing CLI commands to YAML:

1. Identify all `--flag value` pairs
2. Convert to `flag: value` in YAML
3. Lists use YAML list syntax:
   ```yaml
   modules:
     - module1
     - module2
   ```
4. Tuples can be lists or comma-separated strings:
   ```yaml
   initial_procs: [2, 2, 2]
   # OR
   initial_procs: "2,2,2"
   ```

## Next Steps

- See `run.yaml` for a fully documented template
- Check `GETTING_STARTED.md` for general usage
- Read `QUICK_START_SCALING.md` for scaling test basics
- Explore `examples/` for more YAML configurations

## Summary

YAML configuration makes HPC-ScaleTest:
- ✅ **Easier to use** - No more super-long command lines
- ✅ **More maintainable** - Track configurations in version control
- ✅ **More flexible** - Mix YAML and CLI as needed
- ✅ **Fully automated** - Auto-detection handles system-specific details
- ✅ **Portable** - Same config works across different HPC systems

Happy scaling! 🚀
