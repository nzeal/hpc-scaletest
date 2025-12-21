# Hardware Topology Detection System

## Overview

HPC-ScaleTest v2.0.0 includes a centralized, portable hardware topology detection
system that automatically infers hardware configuration and generates correct MPI
launch commands without hardcoded or site-specific constants.

## Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         User Request: Run on N nodes                         │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    SLURM Scheduler (backends/schedulers/slurm.py)           │
│  • Calls topology detector                                                   │
│  • Generates #SBATCH directives                                             │
│  • Includes GPU binding script                                               │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    Topology Detector (core/topology.py)                      │
│  • Detects CPU cores per node                                               │
│  • Detects GPUs per node                                                    │
│  • Computes MPI mapping                                                      │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
              ┌───────────────────────┼───────────────────────┐
              │                       │                       │
              ▼                       ▼                       ▼
┌──────────────────────┐ ┌──────────────────────┐ ┌──────────────────────┐
│   SLURM Env Vars     │ │   sinfo/scontrol     │ │  System Introspection │
│  (runtime detection) │ │  (pre-job detection) │ │  (/proc, nvidia-smi)  │
└──────────────────────┘ └──────────────────────┘ └──────────────────────┘
```

## Detection Methods

### Priority 1: SLURM Environment Variables (Runtime)

When running inside a SLURM job, these environment variables are used:

| Variable | Purpose |
|----------|---------|
| `SLURM_CPUS_ON_NODE` | CPU cores on current node |
| `SLURM_JOB_CPUS_PER_NODE` | CPUs per node (format: "32" or "32(x4)") |
| `SLURM_GPUS_ON_NODE` | GPUs on current node |
| `SLURM_GPUS_PER_NODE` | Requested GPUs per node |
| `SLURM_JOB_PARTITION` | Partition name |

This is the most reliable method because it reflects the actual allocated resources.

### Priority 2: SLURM Queries (Pre-submission)

On login nodes before job submission, partition information is queried:

```bash
# Query partition hardware
sinfo -p PARTITION -N -h -o "%c %m %G"
# Returns: cores memory gres(gpu:type:count)

# Query partition configuration
scontrol show partition PARTITION
```

### Priority 3: System Introspection (Fallback)

Direct system inspection is used as fallback:

| Method | Information |
|--------|-------------|
| `lscpu` | CPU cores, sockets, threads |
| `/proc/cpuinfo` | CPU topology |
| `nvidia-smi` | NVIDIA GPU count and model |
| `rocm-smi` | AMD GPU count and model |

## MPI Mapping Algorithm

### For GPU Jobs

When GPUs are present, the mapping follows this algorithm:

```
ranks_per_node = gpus_per_node     # 1 MPI rank per GPU
cores_per_rank = cpu_cores / gpus  # Distribute cores evenly
gpus_per_rank = 1                  # Each rank gets one GPU
```

**Example: 32 cores, 4 GPUs**
```
ranks_per_node = 4
cores_per_rank = 32 / 4 = 8
gpus_per_rank = 1
```

### For CPU Jobs

When no GPUs are present:

```
ranks_per_node = physical_cores    # 1 MPI rank per core
cores_per_rank = 1
gpus_per_rank = 0
```

## MPI Command Generation

### OpenMPI (Correct Syntax)

**CORRECT:**
```bash
mpirun -np 16 \
    --map-by ppr:4:node \
    --bind-to core \
    --cpus-per-proc 8 \
    ./bind_gpu.sh ./app
```

**INCORRECT (DO NOT USE):**
```bash
# This syntax is NOT valid for OpenMPI
mpirun -np 16 --map-by ppr:4:node:PE=8 ./app
```

The `:PE=Y` suffix is not supported by OpenMPI. The correct way to specify
CPUs per process is via the separate `--cpus-per-proc` option.

### Intel MPI

```bash
mpirun -np 16 -ppn 4 ./bind_gpu.sh ./app
```

Note: Intel MPI does NOT support `--report-bindings` (will cause failure).

### Cray MPICH (srun)

```bash
srun -n 16 --ntasks-per-node 4 --cpus-per-task 8 ./app
```

## GPU Binding

GPU binding is handled via a wrapper script that sets `CUDA_VISIBLE_DEVICES`
based on the local MPI rank:

```bash
#!/bin/bash
# Determine local rank
if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
    LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ -n "$PMI_LOCAL_RANK" ]; then
    LOCAL_RANK=$PMI_LOCAL_RANK
elif [ -n "$SLURM_LOCALID" ]; then
    LOCAL_RANK=$SLURM_LOCALID
fi

# Bind to GPU
export CUDA_VISIBLE_DEVICES=$LOCAL_RANK
exec "$@"
```

Supported environment variables for local rank detection:

| MPI Implementation | Variable |
|-------------------|----------|
| OpenMPI | `OMPI_COMM_WORLD_LOCAL_RANK` |
| Intel MPI | `MPI_LOCALRANKID` |
| MPICH/Hydra | `PMI_LOCAL_RANK` |
| SLURM | `SLURM_LOCALID` |
| MVAPICH2 | `MV2_COMM_WORLD_LOCAL_RANK` |

## Generated SLURM Job Script

For a 4-node GPU job on a system with 32 cores and 4 GPUs per node:

```bash
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --partition=gpu_partition
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:4
#SBATCH --time=02:00:00
#SBATCH --job-name=scaling_test

# Environment setup
module load cuda openmpi

# Set OpenMP threads
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Create GPU binding script
cat > bind_gpu.sh << 'EOF'
#!/bin/bash
if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
    export CUDA_VISIBLE_DEVICES=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ -n "$SLURM_LOCALID" ]; then
    export CUDA_VISIBLE_DEVICES=$SLURM_LOCALID
fi
exec "$@"
EOF
chmod +x bind_gpu.sh

# Run application
mpirun -np 16 \
    --map-by ppr:4:node \
    --bind-to core \
    --cpus-per-proc 8 \
    ./bind_gpu.sh ./application
```

## Documented Assumptions

The following assumptions are made and documented:

1. **GPU jobs use 1 MPI rank per GPU**: This is the standard configuration for
   GPU-accelerated HPC applications. Override with `user_ranks_per_node` if needed.

2. **CPUs are evenly divisible by GPUs**: The mapping algorithm divides
   `cpu_cores / gpus`. Non-integer results are truncated.

3. **Hyperthreading is typically disabled on HPC systems**: The default
   `threads_per_core=1` reflects this. Detection methods will update this
   if hyperthreading is detected.

4. **OpenMPI 4.x+ is assumed for `--cpus-per-proc`**: Earlier versions may
   not support this option. The code checks `supports_cpus_per_proc` before
   using it.

## API Usage

### Basic Detection

```python
from core.topology import detect_topology, compute_mpi_mapping

# Detect topology (auto-selects best method)
topology = detect_topology(partition="gpu_partition")
print(f"Cores: {topology.cpu_cores}, GPUs: {topology.gpus}")

# Compute mapping
topology, mapping = compute_mpi_mapping(num_nodes=4, partition="gpu_partition")
print(f"Ranks/node: {mapping.ranks_per_node}, Cores/rank: {mapping.cores_per_rank}")
```

### Advanced Usage

```python
from core.topology import TopologyDetector, NodeTopology, GPUVendor
from core.mpi_command import MPICommandGenerator, MPIInfo, MPIImplementation

# Create custom topology
topology = NodeTopology(
    cpu_cores=32,
    gpus=4,
    gpu_vendor=GPUVendor.NVIDIA
)

# Compute mapping with user overrides
detector = TopologyDetector()
mapping = detector.compute_mpi_mapping(
    topology=topology,
    num_nodes=4,
    user_ranks_per_node=2,  # Override: 2 ranks instead of 4
    user_cores_per_rank=16, # Override: 16 cores instead of 8
)

# Generate MPI command
generator = MPICommandGenerator()
cmd = generator.generate(
    topology=topology,
    mapping=mapping,
    executable="./app",
    num_nodes=4,
    gpu_binding_script="./bind_gpu.sh"
)
print(" ".join(cmd))
```

## Testing

Run the topology detection test:

```bash
python3 tests/test_topology_detection.py --verbose
```

Expected output for Leonardo Booster-style configuration:
```
Example: 32 CPU cores, 4 A100 GPUs
  Mapping: 4 ranks/node, 8 cores/rank
  Command: mpirun -np 16 --map-by ppr:4:node --bind-to core --cpus-per-proc 8 ./app
  ✓ PASS: Correct mapping computed
```
