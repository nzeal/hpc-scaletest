# HPC Auto - Quick Start Guide

Get started with automated HPC scaling tests in 5 minutes!

---

## Installation

```bash
# Navigate to the framework directory
cd hpc-scaletest2

# No additional dependencies required!
# The framework uses Python standard library + existing dependencies
```

---

## Your First Scaling Test (3 Steps)

### Step 1: Prepare Your Code

Make sure you have either:
- **Local code**: A directory with your HPC application
- **Git repository**: A URL to clone your application

### Step 2: Run the Test

```bash
# For local code
python hpc_auto.py /path/to/your/code --scaling strong --nodes 4

# For Git repository
python hpc_auto.py https://github.com/user/repo.git --scaling strong --nodes 4
```

### Step 3: View Results

```bash
# Reports are generated in the output directory
cat output/*/StrongScalingReport.txt
```

**That's it!** 🎉

---

## Example Workflows

### Example 1: Strong Scaling Test (Fixed Problem Size)

```bash
python hpc_auto.py /home/user/my-hpc-app \
    --scaling strong \
    --nodes 8 \
    --procs-per-node 128 \
    --verbose
```

**Output:**
```
===============================================================
Strong Scaling Efficiency Report
===============================================================
Nodes    Procs     Time(s)    Speedup    Efficiency
---------------------------------------------------------------
1        128       100.00     1.00       100.0%
2        256       52.00      1.92       96.2%
4        512       28.00      3.57       89.3%
8        1024      16.00      6.25       78.1%
===============================================================
```

---

### Example 2: Weak Scaling Test (Growing Problem Size)

```bash
python hpc_auto.py /home/user/my-hpc-app \
    --scaling weak \
    --nodes 16 \
    --initial-domain "10.0,10.0,10.0" \
    --initial-cells "256,256,256" \
    --verbose
```

**What happens:**
- 1 node: 256³ cells, domain 10×10×10
- 2 nodes: 512³ cells (scaled), domain 20×10×10 (scaled)
- 4 nodes: 512×512×256 cells, domain 20×20×10
- ... problem size grows with node count

---

### Example 3: GPU Test

```bash
python hpc_auto.py /home/user/cuda-app \
    --hardware gpu \
    --nodes 4 \
    --procs-per-node 4 \
    --gpus-per-node 1 \
    --partition gpu_prod
```

---

### Example 4: Clone from Git and Test

```bash
python hpc_auto.py https://github.com/user/hpc-app.git \
    --scaling strong \
    --nodes 8 \
    --time-limit "04:00:00" \
    --account "my_project"
```

---

### Example 5: Generate Scripts Only (No Auto-Submit)

```bash
python hpc_auto.py /home/user/my-app \
    --scaling strong \
    --nodes 4 \
    --no-submit
```

Then submit manually:
```bash
cd output/my-app_strong_*/
sbatch nodes_1/job.sh
sbatch nodes_2/job.sh
sbatch nodes_4/job.sh
```

---

## Python API Examples

### Simple API

```python
from engine.orchestrator import create_simple_workflow
from utils.logging_config import setup_logging

# Setup logging
setup_logging()

# Create and run workflow
orchestrator = create_simple_workflow(
    source="/path/to/code",
    scaling_type="strong",
    max_nodes=8,
    hardware_type="cpu"
)

success = orchestrator.run()
print(f"Test completed: {'SUCCESS' if success else 'FAILED'}")
```

### Advanced API with Full Configuration

```python
from engine.orchestrator import HPCOrchestrator, OrchestratorConfig
from pathlib import Path
from utils.logging_config import setup_logging

setup_logging()

# Create detailed configuration
config = OrchestratorConfig(
    source="/path/to/code",
    scaling_type="strong",
    max_nodes=16,
    initial_procs=(4, 4, 2),
    procs_per_node=128,
    time_limit="02:00:00",
    partition="production",
    account="my_project",
    modules=["gcc/11.2.0", "openmpi/4.1.1"],
    auto_submit_jobs=True,
    output_dir=Path("results")
)

# Run workflow
orchestrator = HPCOrchestrator(config)
success = orchestrator.run()
```

---

## Understanding the Workflow

### What Happens Behind the Scenes?

```
1. CODE ACQUISITION
   └─→ Clone Git repo or use local path

2. CODE ANALYSIS
   └─→ Scan README for dependencies
   └─→ Detect build system (CMake, Make, etc.)
   └─→ Identify required modules (MPI, compilers, etc.)

3. COMPILATION
   └─→ Load detected modules
   └─→ Configure build system
   └─→ Compile with appropriate flags

4. TEST GENERATION
   └─→ Generate job configurations
   └─→ Create job scripts (job.sh)
   └─→ Scale input files

5. JOB SUBMISSION
   └─→ Submit to scheduler (or prepare for manual submission)
   └─→ Monitor job execution

6. REPORT GENERATION
   └─→ Parse timing data
   └─→ Calculate speedup & efficiency
   └─→ Generate reports (text + JSON)
```

---

## Output Structure

After running a test, you'll find:

```
output/
└── my_app_strong_20251025_104500/
    ├── nodes_1/
    │   ├── job.sh              # Job submission script
    │   ├── input.dat           # Application input file
    │   ├── metadata.json       # Job metadata
    │   └── out_*.out           # Job output (after execution)
    ├── nodes_2/
    │   └── ...
    ├── nodes_4/
    │   └── ...
    ├── summary.json            # Overall test summary
    ├── StrongScalingReport.txt # Human-readable report
    └── strong_scaling_report.json  # Machine-readable report
```

---

## Common Options

### Scaling Configuration

```bash
--scaling strong      # Fixed problem size, increasing parallelism
--scaling weak        # Growing problem size with parallelism
--nodes 8             # Maximum number of nodes to test
--initial-procs "2,2,2"  # Starting process decomposition (x,y,z)
```

### Hardware Configuration

```bash
--hardware cpu               # CPU-only tests
--hardware gpu               # GPU-enabled tests
--procs-per-node 128         # Processes per node
--gpus-per-node 1            # GPUs per node (for --hardware gpu)
```

### Resource Configuration

```bash
--time-limit "02:00:00"      # Job time limit
--partition "production"     # Partition/queue name
--account "project123"       # Account/project name
```

### Backend Configuration

```bash
--scheduler slurm            # Job scheduler (slurm/pbs/local)
--launcher srun              # MPI launcher (srun/mpirun/mpiexec)
--module-system lmod         # Module system (lmod/tmod/nomod)
```

### Build Configuration

```bash
--build-system cmake         # Build system (cmake/make/autotools)
--modules "gcc/11.2.0,openmpi/4.1.1"  # Modules to load
--no-build                   # Skip build, use existing executable
```

### Behavior Options

```bash
--no-submit       # Generate scripts but don't submit
--no-reports      # Skip report generation
--cleanup         # Clean up Git clones after test
--verbose         # Verbose logging
--debug           # Debug logging
```

---

## Tips for Success

### 1. Improve Auto-Detection with Good README

The better your README, the better auto-detection works:

```markdown
# My HPC App

## Dependencies
- MPI (OpenMPI 4.1+)
- HDF5 1.12+
- CUDA 11.0+ (optional)

## Modules
module load gcc/11.2.0
module load openmpi/4.1.1

## Building
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
```

### 2. Start with Small Tests

```bash
# Test with 2 nodes first
python hpc_auto.py /path/to/code --nodes 2 --verbose

# Then scale up
python hpc_auto.py /path/to/code --nodes 16
```

### 3. Use --no-submit for Verification

```bash
# Generate scripts first
python hpc_auto.py /path/to/code --nodes 4 --no-submit

# Review generated scripts
cat output/*/nodes_1/job.sh

# Submit manually if satisfied
sbatch output/*/nodes_1/job.sh
```

### 4. Specify Modules if Auto-Detection Fails

```bash
python hpc_auto.py /path/to/code \
    --modules "gcc/11.2.0,openmpi/4.1.1,hdf5/1.12.0" \
    --nodes 8
```

### 5. Save Configurations for Reproducibility

```python
# Save your configuration in a Python script
from engine.orchestrator import OrchestratorConfig, HPCOrchestrator

config = OrchestratorConfig(
    source="/path/to/code",
    scaling_type="strong",
    max_nodes=16,
    modules=["gcc/11.2.0", "openmpi/4.1.1"],
    # ... other settings
)

# Version control this file for reproducibility!
orchestrator = HPCOrchestrator(config)
orchestrator.run()
```

---

## Troubleshooting

### Problem: "Build failed"

**Solution 1:** Check if modules are available
```bash
module avail gcc
module avail openmpi
```

**Solution 2:** Manually specify modules
```bash
python hpc_auto.py /path/to/code --modules "gcc/11.2.0,openmpi/4.1.1"
```

**Solution 3:** Use pre-built executable
```bash
python hpc_auto.py /path/to/code --no-build
```

---

### Problem: "No timing data in report"

**Solution:** Ensure your application outputs timing information:
- Pattern: "Time: 123.45"
- Pattern: "Elapsed: 123.45 seconds"
- Pattern: "Wall time: 123.45"

---

### Problem: "Job submission failed"

**Solution 1:** Verify scheduler is available
```bash
which sbatch  # For Slurm
which qsub    # For PBS
```

**Solution 2:** Check partition and account
```bash
sinfo  # List available partitions
sacctmgr show account  # List accounts
```

**Solution 3:** Use local scheduler for testing
```bash
python hpc_auto.py /path/to/code --scheduler local --nodes 1
```

---

## Next Steps

1. **Read Full Documentation**: [Automated Workflow Guide](AUTOMATED_WORKFLOW.md)
2. **Explore Examples**: Check `examples/` directory
3. **Customize Tests**: Learn about OrchestratorConfig options
4. **Analyze Results**: Understand efficiency reports

---

## Quick Reference

### Minimal Command
```bash
python hpc_auto.py /path/to/code --nodes 8
```

### Common Command
```bash
python hpc_auto.py /path/to/code \
    --scaling strong \
    --nodes 16 \
    --partition production \
    --account my_project
```

### Full Command
```bash
python hpc_auto.py /path/to/code \
    --scaling strong \
    --nodes 32 \
    --procs-per-node 128 \
    --initial-procs "4,4,2" \
    --time-limit "04:00:00" \
    --partition production \
    --account my_project \
    --modules "gcc/11.2.0,openmpi/4.1.1,hdf5/1.12.0" \
    --build-system cmake \
    --verbose
```

---

## Getting Help

```bash
# Show all options
python hpc_auto.py --help

# Enable verbose logging
python hpc_auto.py /path/to/code --verbose

# Enable debug logging
python hpc_auto.py /path/to/code --debug
```

---

**Happy scaling! 🚀**
