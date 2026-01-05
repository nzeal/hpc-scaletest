# Changelog

## [2.0.0] - 2024 - Major Refactoring Release

### New: Partition-Aware Hardware Topology Detection

A centralized, partition-aware hardware topology detection system that automatically
infers CPU cores and GPUs per node from SLURM partition configuration, then derives
the correct MPI mapping without hardcoded values.

#### Key Features

- **Automatic Detection**: Queries SLURM partition for hardware topology
- **Derived Mapping**: `ranks_per_node = GPUs`, `cores_per_rank = CPUs/GPUs`
- **Correct Syntax**: Generates `--map-by ppr:R:node:PE=C` mpirun commands
- **Portable**: Works on Leonardo, LUMI, DGX, and any SLURM system

#### Example: Leonardo Booster

For partition `boost_usr_prod` with 32 cores and 4 GPUs:

```
Detected: 32 CPU cores per node, 4 GPUs per node
Derived:  4 ranks per node, 8 cores per rank
```

SLURM directives (full node allocation):
```bash
#SBATCH --ntasks-per-node=32
#SBATCH --gres=gpu:4
```

Generated mpirun commands:
```bash
# 1 node:
mpirun -np 4 --map-by ppr:4:node:PE=8 --report-bindings ./bind.sh $BINARY/iPIC3D os-stdin

# 2 nodes:
mpirun -np 8 --map-by ppr:4:node:PE=8 --report-bindings ./bind.sh $BINARY/iPIC3D os-stdin

# 4 nodes:
mpirun -np 16 --map-by ppr:4:node:PE=8 --report-bindings ./bind.sh $BINARY/iPIC3D os-stdin
```

Formula:
```
-np = (nodes) × (GPUs per node)
--map-by ppr:(GPUs per node):node:PE=(cores per rank)
```

#### Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                    SLURM Job Script Generation                   │
│                  (backends/schedulers/slurm.py)                 │
├─────────────────────────────────────────────────────────────────┤
│                 Partition Topology Detection                     │
│              (core/partition_topology.py)                       │
├─────────────────────────────────────────────────────────────────┤
│  SLURM Env Vars  │  sinfo query  │  scontrol query              │
└─────────────────────────────────────────────────────────────────┘
```

#### Detection Methods (Priority Order)

1. **Slurm Environment Variables** (runtime, inside job)
   - `SLURM_CPUS_ON_NODE`: CPUs on current node
   - `SLURM_GPUS_ON_NODE`: GPUs on current node
   - `SLURM_JOB_PARTITION`: Partition name
   - Most reliable when running inside an allocation

2. **Slurm Queries** (pre-submission, login node)
   - `sinfo -p PARTITION -o "%c %G %m"`: Cores, memory, GRES
   - `scontrol show partition PARTITION`: Detailed config
   - Works before job submission

3. **System Introspection** (fallback)
   - `/proc/cpuinfo`: CPU cores, sockets
   - `nvidia-smi`: NVIDIA GPUs
   - `rocm-smi`: AMD GPUs
   - `lscpu`: Thread count

#### New Core Modules

| Module | Description |
|--------|-------------|
| `core/topology.py` | Hardware topology detection (CPUs, GPUs per node) |
| `core/mpi_command.py` | MPI command generation with correct syntax |

#### Key Features

1. **Automatic Hardware Detection**
   - Detects CPU cores per node from SLURM env vars or `/proc/cpuinfo`
   - Detects GPUs per node via `nvidia-smi`, `rocm-smi`, or SLURM
   - Supports NVIDIA, AMD, and Intel GPUs

2. **Correct MPI Command Generation**
   - OpenMPI: `--map-by ppr:N:node --bind-to core --cpus-per-proc M`
   - Intel MPI: `-ppn N` (no `--report-bindings`)
   - MPICH/Cray: Native srun syntax
   - **FIXED**: The incorrect `ppr:X:node:PE=Y` syntax is no longer used

3. **Automatic GPU Binding**
   - Uses `OMPI_COMM_WORLD_LOCAL_RANK` / `PMI_LOCAL_RANK` / `SLURM_LOCALID`
   - Sets `CUDA_VISIBLE_DEVICES` for NVIDIA GPUs
   - Sets `ROCR_VISIBLE_DEVICES` for AMD GPUs
   - Sets `ZE_AFFINITY_MASK` for Intel GPUs

4. **Mapping Algorithm**
   - GPU jobs: `ranks_per_node = gpus_per_node` (1 rank per GPU)
   - `cores_per_rank = cpu_cores / ranks_per_node`
   - CPU jobs: `ranks_per_node = physical_cores`

#### Example: Leonardo Booster (32 CPUs, 4 GPUs)

The system automatically derives:
```bash
mpirun -np 16 \
    --map-by ppr:4:node \
    --bind-to core \
    --cpus-per-proc 8 \
    ./bind_gpu.sh ./app
```

This works without any site-specific constants - the same code works on any GPU cluster.

### Critical Bug Fixes

#### Fixed: Incorrect OpenMPI `--map-by` Syntax
- **Problem**: The syntax `--map-by ppr:X:node:PE=Y` was incorrect for OpenMPI.
- **Fix**: Changed to use separate `--map-by ppr:X:node --bind-to core --cpus-per-proc Y`.
- **Files Fixed**:
  - `backends/launchers/mpirun.py`
  - `utils/advanced_gpu_manager.py`
  - `utils/hardware_autodetect.py`
  - `tests/test_gpu_scaling.py`

### Removed Duplicate/Versioned Files

The following redundant files were removed to eliminate technical debt:

#### Source Manager Variants
- ❌ `utils/source_manager_HARDENED.py`
- ❌ `utils/source_manager_fixed.py`
- ❌ `utils/source_manager_improved.py`

#### Build Discovery Variants
- ❌ `utils/build_discovery_BUGFIXED.py`
- ❌ `utils/build_discovery_REFINED.py`
- ❌ `utils/build_discovery_fixed.py`
- ❌ `utils/build_discovery_improved.py`

#### Binary Detector Variants
- ❌ `utils/binary_detector_IMPROVED.py`
- ❌ `utils/binary_detector_fixed.py`

#### GPU Scaling Manager Variants
- ❌ `utils/gpu_scaling_manager_HARDENED.py`
- ❌ `utils/gpu_scaling_manager_fixed.py`
- ❌ `utils/gpu_scaling_manager_improved.py`

#### Other Removed Files
- ❌ `utils/advanced_gpu_manager_py37.py.bak`
- ❌ `backends/builds/cmake_BUGFIXED.py`
- ❌ `backends/launchers/mpirun.py.bak`
- ❌ `leonardo_partitions.py` (moved to external config)

#### Removed Tests for Deleted Modules
- ❌ `tests/test_binary_detector_improved.py`
- ❌ `tests/test_fixed_modules.py`
- ❌ `tests/test_improved_modules.py`
- ❌ `tests/test_gpu_scaling_manager_hardened.py`
- ❌ `tests/test_source_manager_hardened.py`
- ❌ `tests/test_build_discovery_refined.py`
- ❌ `tests/test_cuda_bugfix.py`

### Removed Hardcoded Assumptions

#### System-Specific Defaults
- **Before**: `DEFAULT_PROCS_PER_NODE = 128` (Leonardo-specific)
- **After**: Auto-detected from system or configurable via `HPC_SCALETEST_PROCS_PER_NODE`

#### Partition Names
- **Before**: Hardcoded `'X_usr_prod'` default partition
- **After**: Must be specified or set via `HPC_SCALETEST_PARTITION`

#### Application-Specific Parameters
- **Before**: Hardcoded iPIC3D parameter names (`Lx`, `XLEN`, `npcelx`, etc.)
- **After**: Configurable via `variable_map` with support for multiple aliases

### New Features

#### External System Configuration
- Added `config/systems.yaml` for external HPC system definitions
- Added `utils/yaml_system_loader.py` for loading configurations
- Systems can be selected via `HPC_SCALETEST_SYSTEM` environment variable

#### Improved Types Module
- Added `get_default_procs_per_node()` for dynamic CPU detection
- All defaults now configurable via environment variables:
  - `HPC_SCALETEST_PROCS_PER_NODE`
  - `HPC_SCALETEST_TIME_LIMIT`
  - `HPC_SCALETEST_OUTPUT_DIR`
  - `HPC_SCALETEST_LARGE_JOB_THRESHOLD`

#### Enhanced Input File Parser
- Now supports configurable parameter name mappings
- Supports multiple aliases for the same parameter
- Added `find_parameter()` method for easy parameter lookup

### Configuration

#### Environment Variables
| Variable | Description | Default |
|----------|-------------|---------|
| `HPC_SCALETEST_SYSTEM` | System configuration name | auto-detect |
| `HPC_SCALETEST_PARTITION` | SLURM partition | (required) |
| `HPC_SCALETEST_PROCS_PER_NODE` | CPU cores per node | auto-detect |
| `HPC_SCALETEST_TIME_LIMIT` | Default time limit | `01:00:00` |
| `HPC_SCALETEST_OUTPUT_DIR` | Output directory | `output` |
| `HPC_SCALETEST_CONFIG` | Path to systems.yaml | (auto-search) |

### Migration Guide

#### For Users with Custom System Configurations

1. Move system-specific settings to `config/systems.yaml` or `~/.hpc-scaletest/systems.yaml`
2. Set `HPC_SCALETEST_PARTITION` environment variable
3. Update any scripts that relied on deleted versioned modules

#### For Applications Using Non-iPIC3D Parameter Names

Configure custom parameter names in your YAML config:

```yaml
scaling:
  variable_map:
    domain:
      x: [my_domain_x, Lx]
      y: [my_domain_y, Ly]
      z: [my_domain_z, Lz]
    cells:
      x: [my_nx, nxc]
      y: [my_ny, nyc]
      z: [my_nz, nzc]
```

### File Count Reduction

- **Before**: 130 Python files
- **After**: ~110 Python files
- **Removed**: 20 duplicate/obsolete files

### Testing

All remaining tests updated to use correct OpenMPI syntax and work without system-specific assumptions.
