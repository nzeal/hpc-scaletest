# HPC-ScaleTest Fixes Documentation

## Design Philosophy

HPC-ScaleTest is an automated framework designed to generate and execute scaling studies on high-performance computing systems. Its primary goal is to enable reproducible performance evaluation across a wide range of platforms without requiring manual intervention.

### Core Design Principles

1. **SYSTEM AGNOSTIC**: Works on any HPC system (Leonardo, LUMI, Piz Daint, etc.)
2. **COMPILER AGNOSTIC**: No assumptions about compilers or toolchains
3. **APPLICATION AGNOSTIC**: Works with any MPI application
4. **NO HARDCODED VALUES**: All configuration discovered or provided at runtime
5. **FULL AUTOMATION**: From configuration through execution and data collection

### Configuration Sources (Priority Order)

1. **YAML configuration file** (`--config run.yaml`)
2. **Command-line arguments** (`--nodes 16 --partition gpu_prod`)
3. **System auto-detection** (SLURM queries, hardware detection)

**No hardcoded fallback defaults** - if a required value cannot be determined, the framework reports an error with clear guidance.

---

## Required Configuration Fields

### For All Jobs
| Field | Description | Source |
|-------|-------------|--------|
| `max_nodes` | Maximum nodes for scaling tests | YAML or CLI |
| `partition` | SLURM partition name | YAML or CLI |
| `initial_procs` | MPI decomposition [x, y, z] | YAML or CLI |

### For GPU Jobs (hardware_type: gpu)
| Field | Description | Source |
|-------|-------------|--------|
| `gpus_per_node` | GPUs per node | YAML, CLI, or detection |

### Auto-Detected (Optional Override)
| Field | Description | Detection Method |
|-------|-------------|------------------|
| `procs_per_node` | CPU cores per node | SLURM sinfo/scontrol |
| `memory_per_node_gb` | Memory per node | SLURM node info |
| `scheduler` | Job scheduler type | Environment detection |

---

## Fixes Applied

### Fix 1: GPU mpirun Command Format

**Problem**: GPU mpirun command was broken - used incorrect wrapper script format.

**Root Cause**: `slurm.py` line ~454 generated `./gpu_wrapper.sh` instead of `./bind.sh`.

**Solution**: Fixed command format to:
```bash
mpirun -np <total_ranks> --map-by ppr:<ranks_per_node>:node:PE=<cores_per_rank> ./bind.sh $BINARY/<executable> <input>
```

**Example** (4 nodes, 4 GPUs/node, 32 CPUs/node):
```bash
mpirun -np 16 --map-by ppr:4:node:PE=8 ./bind.sh $BINARY/iPIC3D os-stdin
```

**Files Modified**:
- `backends/schedulers/slurm.py`: Fixed GPU command generation
- `backends/launchers/mpirun.py`: Uses bind.sh wrapper

---

### Fix 2: max_nodes YAML Parsing

**Problem**: `max_nodes: 16` in YAML was ignored, defaulting to 4 nodes.

**Root Cause**: `config_parser.py` only recognized `nodes:` key, not `max_nodes:`.

**Solution**: Added support for both keys:
```python
if 'max_nodes' in self.config_data:
    config['max_nodes'] = int(self.config_data['max_nodes'])
elif 'nodes' in self.config_data:
    config['max_nodes'] = int(self.config_data['nodes'])
```

**Files Modified**:
- `utils/config_parser.py`: Added max_nodes + hardware_type + gpus_per_node parsing

---

### Fix 3: GPU Mode Validation

**Problem**: Strong scaling validator used CPU cores instead of GPUs for MPI rank calculation.

**Solution**: For GPU jobs, use gpus_per_node for effective processes per node:
```python
if self.config.hardware_type.lower() == "gpu" and gpus > 0:
    effective_procs_per_node = gpus  # 1 MPI rank per GPU
else:
    effective_procs_per_node = self.config.procs_per_node
```

**Files Modified**:
- `engine/orchestrator.py`: Fixed GPU mode validation

---

### Fix 4: No Hardcoded Defaults

**Problem**: Hardcoded values like `max_nodes=4`, `procs_per_node=128`, `partition='X_usr_prod'` were used as fallback defaults.

**Solution**: Removed all hardcoded defaults for hardware-specific values:

**In `hpc_auto.py` (CLI entry point)**:
```python
# Before
parser.add_argument('--nodes', default=4, ...)
parser.add_argument('--procs-per-node', default=128, ...)
parser.add_argument('--partition', default='X_usr_prod', ...)

# After
parser.add_argument('--nodes', default=None, ...)  # REQUIRED
parser.add_argument('--procs-per-node', default=None, ...)  # Auto-detect or YAML
parser.add_argument('--partition', default=None, ...)  # REQUIRED
```

**In `engine/orchestrator.py`**:
```python
@dataclass
class OrchestratorConfig:
    max_nodes: Optional[int] = None  # MUST be set by user
    partition: Optional[str] = None  # MUST be set by user
    procs_per_node: Optional[int] = None  # Auto-detect or user-set
    gpus_per_node: int = 0  # 0 = CPU mode, user sets for GPU
    
    def __post_init__(self):
        if self.max_nodes is None:
            raise ValueError("max_nodes MUST be specified in YAML")
        if self.partition is None:
            raise ValueError("partition MUST be specified in YAML")
```

**In `utils/config_parser.py` (validation)**:
```python
def validate(self):
    # max_nodes is REQUIRED
    if not has_max_nodes:
        raise ValueError("max_nodes MUST be specified")
    
    # partition is REQUIRED
    if 'partition' not in self.config_data:
        raise ValueError("partition MUST be specified")
    
    # GPU jobs require gpus_per_node
    if hw_type.lower() == 'gpu' and not has_gpus_per_node:
        raise ValueError("gpus_per_node MUST be specified for GPU jobs")
```

**Files Modified**:
- `hpc_auto.py`: Removed hardcoded defaults from argparse
- `engine/orchestrator.py`: Made max_nodes, partition required
- `utils/config_parser.py`: Added validation for required fields
- `utils/node_sequence_scaling.py`: max_nodes now required parameter
- `utils/smart_directional_scaling.py`: max_nodes now required parameter
- `utils/advanced_gpu_manager.py`: Removed hardcoded CPU fallback

---

### Fix 5: Modular Job Generators

**Problem**: CPU and GPU job generation mixed in single file.

**Solution**: Created modular architecture in `engine/job_generators/`:

```
engine/job_generators/
├── __init__.py          # Package exports
├── base.py              # JobGeneratorBase abstract class
├── cpu_job_generator.py # CPUJobGenerator
└── gpu_job_generator.py # GPUJobGenerator with bind.sh
```

**Benefits**:
- Clear separation of concerns
- Easier maintenance
- No code duplication
- Consistent interface

---

## Example YAML Configurations

### Generic GPU Strong Scaling Test
```yaml
# HPC-ScaleTest Configuration - System Agnostic
# All values must be provided - no hardcoded defaults

repository: https://github.com/your-org/your-mpi-app.git

# Hardware Configuration (REQUIRED for GPU jobs)
hardware_type: "gpu"
gpus_per_node: 4           # Discover from your system: sinfo -o "%P %G"
procs_per_node: 32         # Optional - auto-detected from SLURM if not specified

# Scaling Configuration (REQUIRED)
scaling: strong
max_nodes: 16              # Your maximum test scale

# Domain Configuration (application-specific)
initial_domain: [5.12, 5.12, 1]
initial_cells: [128, 128, 1]
initial_procs: [2, 2, 1]

# SLURM Configuration (REQUIRED - specific to YOUR system)
partition: "your_gpu_partition"    # Get from: sinfo -s
account: "your_project_account"    # Get from: sacctmgr show user $USER

# Optional overrides
launcher: mpirun
time_limit: "02:00:00"

# Modules (system-specific - discover with: module avail)
modules:
  - cuda/12.x
  - openmpi/4.x
```

### Generic CPU Strong Scaling Test
```yaml
repository: https://github.com/your-org/your-cpu-app.git

# Hardware Configuration
hardware_type: "cpu"
# procs_per_node: auto-detected from SLURM

# Scaling Configuration (REQUIRED)
scaling: strong
max_nodes: 64

# Domain Configuration
initial_procs: [4, 4, 2]

# SLURM Configuration (REQUIRED)
partition: "your_cpu_partition"
account: "your_project"

launcher: srun
```

### Discovering Your System Configuration

Before running HPC-ScaleTest, discover your system's configuration:

```bash
# Find available partitions
sinfo -s

# Find GPU configuration
sinfo -o "%P %G %c %m" -p <partition_name>

# Find your account
sacctmgr show user $USER

# Find available modules
module avail cuda openmpi intel
```

---

## Verification Tests

Run the test suite to verify all fixes:

```bash
python3 test_fixes.py
```

**Expected Output**:
```
======================================================================
 HPC-ScaleTest Fix Verification Tests
======================================================================

TEST 1: max_nodes YAML Parsing
  ✓ PASSED: max_nodes correctly parsed as 16

TEST 2: Node Sequence Generation
  ✓ PASSED: Sequence correctly includes all powers of 2 up to 16

TEST 3: GPU Job Generator - mpirun Command
  ✓ PASSED: Command format is correct

TEST 4: GPU Job Script Generation
  ✓ PASSED: Job script contains all required elements

TEST 5: CPU Job Generator
  ✓ PASSED: CPU command format is correct

TEST 6: Required Field Validation (No Hardcoded Defaults)
  ✓ Correctly rejected config missing max_nodes
  ✓ Correctly rejected config missing partition
  ✓ Correctly rejected GPU config missing gpus_per_node
  ✓ Complete GPU config validated successfully

======================================================================
 ALL TESTS PASSED ✓
======================================================================
```

---

## Error Messages

If required fields are missing, you'll see clear error messages:

```
======================================================================
CONFIGURATION ERROR: Required values not specified
======================================================================

The following required configuration values are missing:
  • max_nodes/nodes (use --nodes or specify in YAML)
  • partition (use --partition or specify in YAML)
  • gpus_per_node for GPU jobs (use --gpus-per-node or specify in YAML)

HPC-ScaleTest does NOT use hardcoded defaults.
All values must come from:
  1. YAML configuration file (--config)
  2. Command-line arguments
  3. System auto-detection (when available)

Example YAML configuration:
  max_nodes: 16
  partition: boost_usr_prod
  procs_per_node: 32
  gpus_per_node: 4  # For GPU jobs
  initial_procs: [2, 2, 1]
======================================================================
```

---

## Summary of Changes

| File | Change |
|------|--------|
| `backends/schedulers/slurm.py` | Fixed GPU mpirun format, uses bind.sh |
| `backends/launchers/mpirun.py` | Uses bind.sh wrapper for GPU jobs |
| `utils/config_parser.py` | Added max_nodes, validation for required fields |
| `engine/orchestrator.py` | Made max_nodes, partition required; integrated build strategy |
| `hpc_auto.py` | Removed hardcoded defaults from argparse |
| `utils/node_sequence_scaling.py` | max_nodes is now required parameter |
| `utils/smart_directional_scaling.py` | max_nodes is now required parameter |
| `utils/advanced_gpu_manager.py` | Removed hardcoded CPU fallback values |
| `engine/job_generators/` | NEW: Modular CPU/GPU job generators |
| `engine/build_strategy.py` | NEW: Automated build configuration selection |
| `test_fixes.py` | Comprehensive verification tests |

---

## New: Automated Build Strategy Module

### Overview

The new `engine/build_strategy.py` module provides fully automated build configuration:

1. **Application Analysis**: Scans source code to detect requirements
   - CUDA/HIP/SYCL GPU code detection
   - MPI usage detection
   - OpenMP detection
   - Library dependencies (HDF5, FFTW, PETSc)

2. **System Detection**: Queries hardware capabilities
   - CPU core count from SLURM
   - GPU count and vendor (NVIDIA, AMD, Intel)
   - Available compilers

3. **Build Strategy Selection**: Automatically selects build configuration
   - CPU-only, GPU-only, or hybrid builds
   - Appropriate compiler selection
   - CMake flag configuration

### Usage

The build strategy is automatically invoked during the orchestrator's analysis phase:

```python
from engine.build_strategy import auto_configure_build

strategy, app_req, sys_caps = auto_configure_build(
    source_dir=Path("/path/to/source"),
    partition="gpu_partition",
    hardware_type="gpu"  # Optional hint
)

print(f"Build type: {strategy.build_type.value}")
print(f"CMake flags: {strategy.cmake_flags}")
print(f"Modules: {strategy.modules}")
```

### Execution Models

| Model | Description | Detection |
|-------|-------------|-----------|
| `cpu` | CPU-only execution | No GPU code detected |
| `gpu` | GPU-accelerated | CUDA/HIP files found |
| `hybrid` | MPI + GPU | Both MPI and GPU code |

### Build Types

| Type | Use Case |
|------|----------|
| `cpu_release` | Pure CPU MPI applications |
| `gpu_cuda` | NVIDIA GPU applications |
| `gpu_hip` | AMD ROCm applications |
| `hybrid_cuda` | MPI + NVIDIA GPU |
| `hybrid_hip` | MPI + AMD GPU |

---

## Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                    HPC-ScaleTest Workflow                        │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ Step 1: CODE ACQUISITION                                         │
│   • Clone from Git repository                                    │
│   • Or use local source path                                     │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ Step 2: CODE ANALYSIS & BUILD CONFIGURATION                      │
│   • Scan source files for requirements (CUDA, MPI, OpenMP)      │
│   • Detect system hardware (CPUs, GPUs, compilers)              │
│   • Select build strategy (CPU/GPU/hybrid)                       │
│   • Configure CMake flags automatically                          │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ Step 3: BUILD                                                    │
│   • Load required modules                                        │
│   • Run cmake with auto-configured flags                         │
│   • Compile with detected compilers                              │
│   • Auto-detect compiled binary                                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ Step 4: SCALING TEST EXECUTION                                   │
│   • Generate job scripts for each node count                     │
│   • CPU path: mpirun -np N $BINARY input                        │
│   • GPU path: mpirun -np N --map-by ppr:X:node:PE=Y ./bind.sh   │
│   • Submit jobs via SLURM                                        │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ Step 5: RESULTS & REPORTING                                      │
│   • Collect job outputs                                          │
│   • Generate scaling reports                                     │
│   • Performance analysis                                         │
└─────────────────────────────────────────────────────────────────┘
```

---

## Contact

For questions or issues, please open an issue on the GitHub repository.
