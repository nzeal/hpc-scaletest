# HPC-ScaleTest v2.0.0 - Architecture Documentation

## Overview

HPC-ScaleTest is an automated framework for performing strong and weak scaling studies on heterogeneous high-performance computing systems. This document describes the architecture, key modules, and design decisions.

## Key Design Principles

1. **No Hardcoded Hardware Assumptions**: All hardware topology is detected at runtime
2. **Partition-Aware Detection**: Topology is queried from SLURM partition configuration
3. **Separation of Concerns**: Clear module boundaries for different responsibilities
4. **Correct MPI Mapping**: Proper `ppr:X:node:PE=Y` syntax with GPU binding

## Module Structure

```
hpc-scaletest-v2.0.0/
├── core/                     # Core modules (single source of truth)
│   ├── hardware.py          # Hardware topology detection
│   ├── slurm_script.py      # SLURM job script generation  
│   ├── mpi_command.py       # MPI command generation
│   ├── config.py            # Configuration dataclasses
│   ├── types.py             # Type definitions and enums
│   └── abstracts.py         # Abstract base classes
├── backends/
│   ├── schedulers/
│   │   ├── slurm.py         # SLURM scheduler backend
│   │   └── local.py         # Local execution backend
│   └── launchers/
│       ├── mpirun.py        # MPI launcher backend
│       └── srun.py          # SRUN launcher backend
├── engine/
│   ├── scaling.py           # Scaling logic (strong/weak)
│   ├── runner.py            # Test runner
│   └── orchestrator.py      # Job orchestration
└── utils/                   # Utility modules
```

## Critical Distinction: SLURM Allocation vs MPI Execution

This is the most important architectural decision in HPC-ScaleTest:

### SLURM Resource Allocation (Full Node)
```bash
#SBATCH --ntasks-per-node=32   # Request ALL CPUs on node
#SBATCH --gres=gpu:4           # Request ALL GPUs on node
```

This ensures:
- Exclusive node access
- Proper resource accounting
- No interference from other jobs

### MPI Execution (Actual Tasks)
```bash
mpirun -np 16 --map-by ppr:4:node:PE=8 ./bind.sh <executable>
```

This runs:
- 4 MPI ranks per node (1 per GPU)
- 8 CPU cores per rank
- GPU binding via CUDA_VISIBLE_DEVICES

### Why This Matters

The previous implementation incorrectly used:
```bash
#SBATCH --ntasks-per-node=4    # WRONG: only requests 4 CPUs
```

This caused resource allocation issues and potential conflicts with other jobs.

## Hardware Topology Detection

### Module: `core/hardware.py`

The `HardwareTopology` dataclass stores detected hardware:

```python
@dataclass
class HardwareTopology:
    partition: str
    cpu_cores_per_node: int
    gpus_per_node: int
    
    # Derived (computed automatically)
    mpi_ranks_per_node: int   # = gpus_per_node for GPU jobs
    cores_per_rank: int       # = cpu_cores / gpus
```

### Detection Priority

1. **SLURM Environment** (inside job):
   - `SLURM_CPUS_ON_NODE`
   - `SLURM_GPUS_ON_NODE`
   - `CUDA_VISIBLE_DEVICES`

2. **sinfo Query** (login node):
   ```bash
   sinfo -p PARTITION -N -h -o "%c %G %m"
   ```

3. **scontrol Query** (fallback):
   ```bash
   scontrol show partition PARTITION
   scontrol show node NODENAME
   ```

4. **System Introspection** (last resort):
   - `/proc/cpuinfo`
   - `nvidia-smi`

## MPI Command Generation

### Module: `core/mpi_command.py`

Uses the `ppr:X:node:PE=Y` syntax:

```python
cmd = [
    'mpirun',
    '-np', str(total_ranks),
    '--map-by', f'ppr:{ranks_per_node}:node:PE={cores_per_rank}',
    '--report-bindings',
]
```

### GPU Binding

The GPU binding script uses `OMPI_COMM_WORLD_LOCAL_RANK`:

```bash
#!/bin/bash
if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
    LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ -n "$SLURM_LOCALID" ]; then
    LOCAL_RANK=$SLURM_LOCALID
else
    LOCAL_RANK=0
fi
export CUDA_VISIBLE_DEVICES=$LOCAL_RANK
exec "$@"
```

## Example: Leonardo Booster

### Hardware
- Partition: `boost_usr_prod`
- CPUs per node: 32
- GPUs per node: 4

### Derived Configuration
```python
mpi_ranks_per_node = 4      # 1 rank per GPU
cores_per_rank = 8          # 32 / 4 = 8
```

### Generated Commands

| Nodes | SLURM Directive | MPI Command |
|-------|-----------------|-------------|
| 1 | `--ntasks-per-node=32` | `mpirun -np 4 --map-by ppr:4:node:PE=8` |
| 2 | `--ntasks-per-node=32` | `mpirun -np 8 --map-by ppr:4:node:PE=8` |
| 4 | `--ntasks-per-node=32` | `mpirun -np 16 --map-by ppr:4:node:PE=8` |

## Scaling Modes

### Weak Scaling

Problem size grows proportionally with resources:
- Domain and grid scale with node count
- MPI decomposition adjusts accordingly
- Work per process remains constant

### Strong Scaling

Problem size remains constant:
- Domain and grid are fixed
- Only MPI decomposition changes
- Work per process decreases

## Portability

The framework works on any system without code changes:

### LUMI-G (128 cores, 8 GPUs)
```
mpi_ranks_per_node = 8
cores_per_rank = 16
mpirun -np 32 --map-by ppr:8:node:PE=16  # 4 nodes
```

### DGX (128 cores, 8 GPUs)
```
mpi_ranks_per_node = 8
cores_per_rank = 16
mpirun -np 64 --map-by ppr:8:node:PE=16  # 8 nodes
```

### CPU-Only (64 cores, 0 GPUs)
```
mpi_ranks_per_node = 64
cores_per_rank = 1
mpirun -np 256 --map-by ppr:64:node:PE=1  # 4 nodes
```

## Remaining Assumptions

1. **GPU jobs use 1 MPI rank per GPU**: This is the standard configuration for GPU-accelerated HPC codes. Can be overridden if needed.

2. **CPUs evenly divisible by GPUs**: Most HPC systems have balanced ratios. Integer division truncates remainder.

3. **OpenMPI 4.x+ syntax**: Uses `ppr:X:node:PE=Y` format. Other MPI implementations may need different syntax.

4. **NVIDIA GPU binding**: Uses `CUDA_VISIBLE_DEVICES`. AMD GPUs need `ROCR_VISIBLE_DEVICES`.

## Testing

Run the comprehensive test suite:

```bash
python tests/test_comprehensive_topology.py
```

This verifies:
- Hardware topology detection
- MPI configuration derivation
- SLURM job script generation
- Correct separation of allocation vs execution
- Portability across different systems

## Changelog Summary

### v2.0.0 Fixes

1. **BUG FIX**: `--ntasks-per-node` now uses full CPU count (32), not MPI ranks (4)
2. **BUG FIX**: MPI mapping uses `ppr:X:node:PE=Y` syntax consistently
3. **BUG FIX**: `BINARY` path is now absolute build directory (not `${PWD}`)
4. **BUG FIX**: Executable is actual binary (not bind.sh)
5. **BUG FIX**: Node sequences dynamically generated (not hardcoded `[1, 2, 4]`)
6. **REFACTOR**: Unified hardware detection in `core/hardware.py`
7. **REFACTOR**: Clean SLURM script generation in `core/slurm_script.py`
8. **REFACTOR**: Incremental result writing in `core/result_writer.py`
9. **REFACTOR**: Unified execution module in `core/unified_execution.py` (single source of truth)
10. **REFACTOR**: Separate CPU (`engine/cpu_execution.py`) and GPU (`engine/gpu_execution.py`) execution
11. **FEATURE**: Automatic partition-aware topology detection
12. **FEATURE**: Portable across Leonardo, LUMI, DGX, CPU-only systems
13. **FEATURE**: Results written immediately to `run_data/` directory
14. **FEATURE**: CSV output for post-processing

### New Modules

| Module | Purpose |
|--------|---------|
| `core/unified_execution.py` | **SINGLE SOURCE OF TRUTH** - Use this for all execution logic |
| `core/result_writer.py` | Incremental result writing to `run_data/` |
| `engine/gpu_execution.py` | GPU job execution with correct binding |
| `engine/cpu_execution.py` | CPU-only job execution |

### Test Suites

```bash
# Comprehensive bug fix verification
python tests/test_comprehensive_bugfix.py

# Topology and command generation
python tests/test_comprehensive_topology.py

# Refactoring requirements
python tests/test_refactoring_requirements.py
```
