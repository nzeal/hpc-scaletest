# HPC Auto - Automated Scaling Framework

## Overview

HPC Auto provides a fully automated end-to-end workflow for HPC application scaling tests. With minimal input, the framework handles code acquisition, dependency detection, compilation, test execution, and performance reporting.

## Key Features

- **🚀 Minimal User Input**: Provide just a code path or Git URL
- **🧠 Intelligent Analysis**: Automatic dependency detection from README files
- **⚙️ Auto-Compilation**: Handles build systems automatically (CMake, Make, Autotools)
- **📊 Scaling Tests**: Automated strong and weak scaling configurations
- **💻 CPU & GPU Support**: Seamless handling of heterogeneous workloads
- **📈 Efficiency Reports**: Beautiful, publication-ready performance reports

---

## Quick Start

### 1. Command-Line Interface (Simplest)

```bash
# Strong scaling test on local code
python hpc_auto.py /path/to/code --scaling strong --nodes 8

# Weak scaling from Git repository
python hpc_auto.py https://github.com/user/hpc-app.git --scaling weak --nodes 16

# GPU-enabled test
python hpc_auto.py /path/to/cuda-app --hardware gpu --nodes 4
```

### 2. Python API (More Control)

```python
from engine.orchestrator import create_simple_workflow

# Create and run workflow
orchestrator = create_simple_workflow(
    source="/path/to/code",
    scaling_type="strong",
    max_nodes=8,
    hardware_type="cpu"
)

success = orchestrator.run()
```

### 3. Full Python Configuration (Maximum Control)

```python
from engine.orchestrator import HPCOrchestrator, OrchestratorConfig
from pathlib import Path

config = OrchestratorConfig(
    source="/path/to/code",
    scaling_type="strong",
    max_nodes=32,
    procs_per_node=128,
    modules=["gcc/11.2.0", "openmpi/4.1.1"],
    auto_submit_jobs=True
)

orchestrator = HPCOrchestrator(config)
orchestrator.run()
```

---

## Workflow Steps

The framework executes the following steps automatically:

### Step 1: Code Acquisition 📥

**What it does:**
- Accepts local path or Git repository URL
- Clones Git repositories to workspace directory
- Validates source code availability

**Examples:**
```bash
# Local path
python hpc_auto.py /home/user/my-hpc-app --scaling strong --nodes 4

# Git HTTPS
python hpc_auto.py https://github.com/user/hpc-app.git --scaling strong --nodes 4

# Git SSH
python hpc_auto.py git@github.com:user/hpc-app.git --scaling strong --nodes 4
```

---

### Step 2: Code Analysis & Dependency Detection 🔍

**What it does:**
- Scans README.md, INSTALL, and other documentation files
- Detects build system (CMake, Make, Autotools, Meson)
- Identifies required modules (MPI, compilers, libraries)
- Detects hardware requirements (CUDA, OpenMP)
- Extracts build commands and flags

**Detection Capabilities:**

| Feature | Detection Method |
|---------|-----------------|
| **Build System** | CMakeLists.txt, Makefile, configure scripts |
| **MPI Required** | Keywords: mpi, mpicc, openmpi, mpich |
| **CUDA Required** | Keywords: cuda, nvcc, cudart |
| **OpenMP Required** | Keywords: openmp, omp, #pragma omp |
| **Compilers** | gcc, intel, clang references |
| **Modules** | Pattern matching: gcc/11.2.0, openmpi/4.1.1 |

**Confidence Scoring:**
The analyzer assigns a confidence score (0-1) based on:
- Build system files found: 30%
- Modules detected: 20%
- Compilers identified: 15%
- Build commands found: 15%
- Dependencies detected: 20%

---

### Step 3: Environment Setup & Compilation ⚙️

**What it does:**
- Loads detected modules automatically
- Sets up compiler environment (MPI, CUDA)
- Runs build system with appropriate flags
- Installs executable to test directory

**Supported Build Systems:**

#### CMake (Default)
```bash
# Auto-detected and executed:
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=mpicc
make -j8
make install
```

#### Make
```bash
# Auto-detected and executed:
make -j8
make install
```

#### Autotools
```bash
# Auto-detected and executed:
./configure --prefix=/install/path
make -j8
make install
```

**Build Flags (Auto-configured):**
- `CMAKE_C_COMPILER=mpicc` (if MPI detected)
- `CMAKE_CXX_COMPILER=mpicxx` (if MPI detected)
- `CMAKE_CUDA_COMPILER=nvcc` (if CUDA detected)
- `-fopenmp` flags (if OpenMP detected)

---

### Step 4: Scaling Test Automation 🧪

**What it does:**
- Generates job configurations for strong/weak scaling
- Creates job submission scripts (job.sh)
- Configures input files with scaled parameters
- Submits jobs to scheduler (or prepares for manual submission)
- Monitors job execution

#### Strong Scaling
**Problem size stays constant, parallelism increases**

```python
# Configuration
scaling_type="strong"
max_nodes=16
initial_procs=(2, 2, 2)  # Start with 8 processes

# Generated configurations:
# - 1 node:  128 procs (8×8×2)
# - 2 nodes: 256 procs (16×8×2)
# - 4 nodes: 512 procs (16×16×2)
# - 8 nodes: 1024 procs (32×16×2)
# - 16 nodes: 2048 procs (32×32×2)
```

#### Weak Scaling
**Problem size grows proportionally with parallelism**

```python
# Configuration
scaling_type="weak"
max_nodes=16
initial_procs=(2, 2, 2)
initial_domain=(10.0, 10.0, 10.0)
initial_cells=(256, 256, 256)

# Generated configurations:
# - 1 node:  8 procs, domain: 10×10×10, cells: 256×256×256
# - 2 nodes: 16 procs, domain: 20×10×10, cells: 512×256×256
# - 4 nodes: 32 procs, domain: 20×20×10, cells: 512×512×256
# (Domain and cells scale with processor count)
```

**Job Script Generation:**

Each configuration generates:
```
output/
└── test_name_strong_20251025_104500/
    ├── nodes_1/
    │   ├── job.sh            # Submission script
    │   ├── input.dat         # Scaled input file
    │   ├── metadata.json     # Job metadata
    │   └── out_*.out         # Job output (after run)
    ├── nodes_2/
    │   └── ...
    └── summary.json          # Overall summary
```

**Example job.sh (Slurm):**
```bash
#!/bin/bash
#SBATCH --job-name=test_nodes_4
#SBATCH --nodes=4
#SBATCH --ntasks=512
#SBATCH --time=02:00:00
#SBATCH --partition=X_usr_prod
#SBATCH --account=cin_X

module load gcc/11.2.0
module load openmpi/4.1.1
export OMP_NUM_THREADS=1

srun -n 512 ./my_app input.dat
```

---

### Step 5: Job Submission & Monitoring 📤

**Automatic Submission (Default):**
```bash
python hpc_auto.py /path/to/code --scaling strong --nodes 8
# Jobs submitted automatically, framework monitors completion
```

**Manual Submission (Prepare Scripts Only):**
```bash
python hpc_auto.py /path/to/code --scaling strong --nodes 8 --no-submit
# Job scripts generated, submit later:
# sbatch output/test_strong_*/nodes_1/job.sh
# sbatch output/test_strong_*/nodes_2/job.sh
# ...
```

**Monitoring:**
- Framework polls job status every 10 seconds
- Handles job states: PENDING, RUNNING, COMPLETED, FAILED
- Collects output files and timing data

---

### Step 6: Post-Processing & Report Generation 📊

**What it does:**
- Parses job outputs for timing data
- Calculates speedup and efficiency metrics
- Generates text and JSON reports

#### Example Strong Scaling Report

```
================================================================================
Strong Scaling Efficiency Report
================================================================================
Test Name: my_hpc_app
Generated: 2025-10-25 10:45:00
================================================================================

Nodes      Procs      Time(s)      Speedup      Efficiency   Decomposition
--------------------------------------------------------------------------------
1          128        100.00       1.00         100.0%       (8×8×2)
2          256        52.00        1.92         96.2%        (16×8×2)
4          512        28.00        3.57         89.3%        (16×16×2)
8          1024       16.00        6.25         78.1%        (32×16×2)
16         2048       10.00        10.00        62.5%        (32×32×2)
================================================================================

Summary Statistics:
----------------------------------------
  Average Efficiency: 85.2%
  Maximum Efficiency: 100.0%
  Minimum Efficiency: 62.5%
  Maximum Speedup: 10.00x
  Parallel Efficiency (max config): 62.5%
================================================================================

Interpretation Guide:
----------------------------------------
  Strong Scaling: Problem size is fixed, parallelism increases
  - Speedup: How much faster compared to baseline configuration
  - Efficiency: How well resources are utilized (100% = ideal)
  - Good strong scaling: Efficiency > 70% up to high node counts
================================================================================
```

#### Report Files Generated

```
output/test_strong_20251025_104500/
├── StrongScalingReport.txt        # Human-readable report
├── strong_scaling_report.json    # Machine-readable data
└── summary.json                   # Raw test data
```

#### JSON Report Structure

```json
{
  "test_name": "my_hpc_app",
  "scaling_type": "strong",
  "generated_at": "2025-10-25T10:45:00",
  "results": [
    {
      "nodes": 1,
      "procs": 128,
      "time_seconds": 100.0,
      "speedup": 1.0,
      "efficiency_percent": 100.0,
      "procs_decomposition": [8, 8, 2]
    },
    ...
  ],
  "statistics": {
    "average_efficiency": 85.2,
    "max_efficiency": 100.0,
    "min_efficiency": 62.5,
    "max_speedup": 10.0
  }
}
```

---

## Usage Examples

### Example 1: Minimal Configuration

```bash
python hpc_auto.py /path/to/code --scaling strong --nodes 8
```

**What happens:**
1. Uses local code at `/path/to/code`
2. Detects build system and dependencies automatically
3. Compiles code with detected configuration
4. Runs strong scaling test up to 8 nodes
5. Generates efficiency report

### Example 2: Git Repository with Custom Resources

```bash
python hpc_auto.py https://github.com/user/app.git \
    --scaling weak \
    --nodes 16 \
    --procs-per-node 64 \
    --partition gpu_prod \
    --account project123 \
    --time-limit 04:00:00
```

### Example 3: GPU Test

```bash
python hpc_auto.py /path/to/cuda-app \
    --hardware gpu \
    --nodes 4 \
    --procs-per-node 4 \
    --gpus-per-node 1 \
    --scaling strong
```

### Example 4: Manual Module Specification

```bash
python hpc_auto.py /path/to/code \
    --scaling strong \
    --nodes 8 \
    --modules "gcc/11.2.0,openmpi/4.1.1,hdf5/1.12.0" \
    --build-system cmake
```

### Example 5: Generate Scripts Without Submission

```bash
python hpc_auto.py /path/to/code \
    --scaling strong \
    --nodes 8 \
    --no-submit
```

Then submit manually:
```bash
cd output/test_strong_*/
sbatch nodes_1/job.sh
sbatch nodes_2/job.sh
# ... etc
```

---

## Configuration Options

### Full CLI Options

```bash
python hpc_auto.py --help
```

| Option | Description | Default |
|--------|-------------|---------|
| `source` | Path or Git URL (required) | - |
| `--scaling` | Scaling type: strong/weak | strong |
| `--nodes` | Maximum nodes to test | 4 |
| `--hardware` | Hardware: cpu/gpu | cpu |
| `--procs-per-node` | Processes per node | 128 |
| `--gpus-per-node` | GPUs per node | 0 |
| `--scheduler` | Scheduler: slurm/pbs/local | slurm |
| `--launcher` | Launcher: srun/mpirun/mpiexec | auto |
| `--module-system` | Module system: lmod/tmod/nomod | lmod |
| `--build-system` | Build system: cmake/make/autotools | auto |
| `--modules` | Comma-separated module list | auto |
| `--partition` | Partition/queue name | X_usr_prod |
| `--account` | Account/project name | cin_X |
| `--time-limit` | Job time limit | 02:00:00 |
| `--no-submit` | Generate scripts only | False |
| `--no-reports` | Skip report generation | False |
| `--verbose` | Verbose logging | False |
| `--debug` | Debug logging | False |

---

## Python API

### Simple API

```python
from engine.orchestrator import create_simple_workflow

orchestrator = create_simple_workflow(
    source="/path/to/code",
    scaling_type="strong",
    max_nodes=8,
    hardware_type="cpu",
    procs_per_node=128,
    auto_submit_jobs=True
)

success = orchestrator.run()
```

### Advanced API

```python
from engine.orchestrator import HPCOrchestrator, OrchestratorConfig
from pathlib import Path

config = OrchestratorConfig(
    source="/path/to/code",
    scaling_type="strong",
    hardware_type="cpu",
    max_nodes=32,
    initial_procs=(4, 4, 2),
    initial_domain=(40.96, 20.48, 1.0),
    initial_cells=(896, 512, 1),
    procs_per_node=128,
    gpus_per_node=0,
    time_limit="04:00:00",
    partition="production",
    account="my_project",
    scheduler="slurm",
    launcher="srun",
    module_system="lmod",
    build_system="cmake",
    modules=["gcc/11.2.0", "openmpi/4.1.1"],
    build_flags={
        "CMAKE_BUILD_TYPE": "Release",
        "ENABLE_HDF5": "ON"
    },
    auto_submit_jobs=True,
    generate_reports=True,
    output_dir=Path("results"),
    test_name="my_test"
)

orchestrator = HPCOrchestrator(config)
success = orchestrator.run()
```

---

## Best Practices

### 1. README Documentation

For best auto-detection results, ensure your README includes:

```markdown
# My HPC Application

## Dependencies
- MPI (OpenMPI 4.1+ or MPICH 3.4+)
- HDF5 1.12+
- CUDA 11.0+ (for GPU version)

## Modules
module load gcc/11.2.0
module load openmpi/4.1.1
module load hdf5/1.12.0

## Building
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
make install
```

### 2. Start Small

Begin with small node counts to verify the workflow:
```bash
python hpc_auto.py /path/to/code --nodes 2 --verbose
```

### 3. Use --no-submit for Testing

Generate scripts first to verify configuration:
```bash
python hpc_auto.py /path/to/code --nodes 4 --no-submit
# Review generated scripts, then submit manually
```

### 4. Weak Scaling Parameters

For weak scaling, provide initial problem size:
```bash
python hpc_auto.py /path/to/code \
    --scaling weak \
    --initial-domain "10.0,10.0,10.0" \
    --initial-cells "256,256,256"
```

### 5. GPU Testing

Adjust processes per node for GPU nodes:
```bash
python hpc_auto.py /path/to/cuda-app \
    --hardware gpu \
    --procs-per-node 4 \
    --gpus-per-node 1
```

---

## Troubleshooting

### Build Fails

**Problem:** Auto-build doesn't work

**Solutions:**
1. Check README documentation quality
2. Manually specify modules:
   ```bash
   --modules "gcc/11.2.0,openmpi/4.1.1"
   ```
3. Manually specify build system:
   ```bash
   --build-system cmake
   ```
4. Use `--no-build` and provide pre-built executable

### Module Loading Issues

**Problem:** Modules not found or not loading

**Solutions:**
1. Verify module names with `module avail`
2. Use exact module names:
   ```bash
   --modules "gcc/11.2.0,openmpi/4.1.1"
   ```
3. Switch module system:
   ```bash
   --module-system tmod
   ```

### Job Submission Fails

**Problem:** Jobs fail to submit

**Solutions:**
1. Check scheduler availability
2. Verify partition and account names
3. Use `--no-submit` to check generated scripts
4. Switch to local scheduler for testing:
   ```bash
   --scheduler local
   ```

### No Timing Data

**Problem:** Reports show no timing information

**Solution:** Ensure your application outputs timing data in a recognizable format. The framework looks for patterns like:
- "Time: 123.45"
- "Elapsed time: 123.45s"
- "Wall time: 123.45"

---

## Integration with Existing Framework

The automated workflow builds on top of the existing HPC-ScaleTest framework:

```
Automated Layer (New)
  ├── hpc_auto.py (CLI)
  ├── orchestrator.py (High-level workflow)
  ├── code_acquisition.py (Git/local)
  ├── readme_analyzer.py (Dependency detection)
  └── report_generator.py (Efficiency reports)

Existing Framework (Unchanged)
  ├── Test definition
  ├── Backend factory
  ├── Schedulers (Slurm, PBS, Local)
  ├── Launchers (srun, mpirun)
  ├── Module systems (Lmod, Tmod)
  └── Build systems (CMake, Make)
```

You can still use the original API for fine-grained control!

---

## Summary

**HPC Auto provides:**
✅ Single command execution  
✅ Intelligent dependency detection  
✅ Automatic compilation  
✅ Strong/weak scaling tests  
✅ CPU and GPU support  
✅ Beautiful efficiency reports  
✅ Flexible configuration  

**Minimal effort, maximum insight!**
