# Quick Start: Running Scaling Tests

## Running Your iPIC3D Strong/Weak Scaling Test

### Option 1: With System Configuration (Recommended)

If you have `leonardo_system.py` in your directory:

```bash
# Strong scaling test (default) - 16 nodes
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git \
    --nodes 16 \
    --partition-name booster \
    --account YOUR_ACCOUNT

# Weak scaling test - 32 nodes
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git \
    --scaling weak \
    --nodes 32 \
    --partition-name booster \
    --account YOUR_ACCOUNT

# With custom launcher (for specific MPI binding)
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git \
    --scaling strong \
    --nodes 16 \
    --partition-name booster \
    --use-custom-launcher mpirun-mapby

# With profiling enabled
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git \
    --nodes 8 \
    --partition-name booster \
    --use-custom-launcher mpirun-nsys
```

**What happens automatically:**
- ✅ System detects 32 CPUs per node from `leonardo_system.py`
- ✅ Partition settings (account, scheduler) auto-loaded
- ✅ Custom launcher (if specified) used for MPI binding
- ✅ Repository cloned and built automatically
- ✅ Scaling jobs generated with correct decomposition
- ✅ Jobs submitted to scheduler

### Option 2: Without System Configuration

If you don't have a system configuration file:

```bash
# Strong scaling test - manually specify everything
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git \
    --scaling strong \
    --nodes 16 \
    --procs-per-node 32 \
    --partition boost_usr_prod \
    --account YOUR_ACCOUNT \
    --time-limit 02:00:00

# Weak scaling test
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git \
    --scaling weak \
    --nodes 32 \
    --procs-per-node 32 \
    --partition boost_usr_prod \
    --account YOUR_ACCOUNT
```

### Option 3: Auto-Detection (Easiest)

If `leonardo_system.py` is in the current directory or parent directory:

```bash
# Just specify the repo and nodes - everything else auto-detected!
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git --nodes 16

# Weak scaling
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git \
    --scaling weak --nodes 32
```

## Understanding the Command Options

### Basic Options

```bash
--scaling strong|weak    # Scaling type (default: strong)
--nodes N                # Maximum number of nodes to test
--partition-name NAME    # Partition from system config (e.g., booster, dcgp)
```

### System Configuration Options

```bash
--system-config FILE     # Path to system config file (auto-detected if not specified)
--partition-name NAME    # Which partition to use from config
--use-custom-launcher L  # Custom launcher name (e.g., mpirun-mapby, mpirun-nsys)
--no-system-config       # Disable system config, use manual settings
```

### Hardware Options

```bash
--hardware cpu|gpu       # Target hardware (default: cpu)
--procs-per-node N       # Override procs per node (auto from system config)
--gpus-per-node N        # GPUs per node (auto from system config for GPU hardware)
```

### Advanced Options

```bash
--initial-procs X,Y,Z    # Initial decomposition (default: 2,2,2)
--initial-domain X,Y,Z   # Domain size for weak scaling
--initial-cells X,Y,Z    # Cell count for weak scaling
--time-limit HH:MM:SS    # Job time limit
--no-submit              # Generate scripts but don't submit
--verbose                # Show detailed progress
```

## What Gets Generated

After running the command, you'll see:

```
output/
  iPIC3D_strong_20241025_133000/
    nodes_1/
      job.sh              # Job script for 1 node
      input.inp           # Modified input file
    nodes_2/
      job.sh
      input.inp
    nodes_4/
      ...
    summary.txt           # Test configuration summary
    efficiency_report.txt # Scaling efficiency results
```

## Checking Job Status

After submission:

```bash
# Check job queue (Slurm)
squeue -u $USER

# View job output
tail -f output/iPIC3D_strong_*/nodes_8/job.out

# Check all jobs
for dir in output/iPIC3D_strong_*/nodes_*/; do
    echo "=== $dir ==="
    tail -20 "$dir/job.out"
done
```

## Example Workflow: Complete Test Run

### 1. Strong Scaling Test with Leonardo Config

```bash
# Run strong scaling from 1 to 16 nodes
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git \
    --scaling strong \
    --nodes 16 \
    --partition-name booster \
    --time-limit 01:00:00 \
    --verbose
```

**Expected output:**
```
HPC Auto - Automated HPC Scaling Framework
============================================================
Attempting to auto-detect system configuration...
✓ Auto-detected system: leonardo
Applying system config from partition: booster
  procs_per_node: 32 (from system config)
  scheduler: slurm (from system config)
  launcher: srun (from system config)
  module_system: tmod4 (from system config)
  partition: boost_usr_prod (from system config)

Cloning repository...
✓ Repository cloned to workspace/iPIC3D-CPU-NS

Building application...
✓ Build completed

Generating scaling configurations...
✓ Generated 5 job configurations:
  - nodes_1: 32 procs, decomposition (4,4,2)
  - nodes_2: 64 procs, decomposition (8,4,2)
  - nodes_4: 128 procs, decomposition (8,8,2)
  - nodes_8: 256 procs, decomposition (16,8,2)
  - nodes_16: 512 procs, decomposition (16,16,2)

Submitting jobs...
✓ Job 12345 submitted: nodes_1
✓ Job 12346 submitted: nodes_2
✓ Job 12347 submitted: nodes_4
✓ Job 12348 submitted: nodes_8
✓ Job 12349 submitted: nodes_16

Results will be saved to: output/iPIC3D_strong_20241025_133000
```

### 2. Weak Scaling Test

```bash
# Run weak scaling - problem size grows with nodes
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git \
    --scaling weak \
    --nodes 32 \
    --partition-name booster \
    --initial-domain 10.0,10.0,10.0 \
    --initial-cells 256,256,256 \
    --verbose
```

### 3. With Custom MPI Binding

```bash
# Use custom launcher for optimal process placement
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git \
    --nodes 16 \
    --partition-name booster \
    --use-custom-launcher mpirun-mapby \
    --verbose
```

This will use the `mpirun-mapby` launcher defined in `leonardo_system.py` with:
```
mpirun -np <N> --map-by socket:PE=8 --rank-by core --report-bindings
```

### 4. Profiling Run

```bash
# Run with NVIDIA Nsight profiling
python hpc_auto.py https://github.com/thenitinshukla/iPIC3D-CPU-NS.git \
    --nodes 4 \
    --partition-name booster \
    --use-custom-launcher mpirun-nsys \
    --verbose
```

This generates per-rank profiling data using:
```
mpirun -np <N> nsys profile -t mpi,cuda --force-overwrite=true \
    -o report_rank%q{OMPI_COMM_WORLD_RANK} --stats=true
```

## Troubleshooting

### If system config is not detected:

```bash
# Explicitly specify the config file
python hpc_auto.py <repo-url> \
    --system-config leonardo_system.py \
    --partition-name booster \
    --nodes 16
```

### If you get "wrong number of procs per node":

Check your system config has correct values:
```python
# In leonardo_system.py
"processor": {
    "num_cpus": 32,  # Make sure this matches your hardware
}
```

### If custom launcher not found:

Make sure the launcher is registered in your system config:
```python
@register_launcher('mpirun-mapby')
class MpirunMapbyLauncher(JobLauncher):
    # ...
```

### To run without submitting (dry-run):

```bash
python hpc_auto.py <repo-url> --nodes 8 --no-submit
# Review generated scripts in output/ directory
# Then submit manually later
```

## Quick Command Reference

```bash
# Simplest: auto-detect everything
python hpc_auto.py <repo-url> --nodes 16

# Strong scaling with system config
python hpc_auto.py <repo-url> --nodes 16 --partition-name booster

# Weak scaling
python hpc_auto.py <repo-url> --scaling weak --nodes 32 --partition-name booster

# Custom launcher
python hpc_auto.py <repo-url> --nodes 8 --use-custom-launcher mpirun-mapby

# GPU test
python hpc_auto.py <repo-url> --hardware gpu --nodes 4 --partition-name booster

# Dry run (don't submit)
python hpc_auto.py <repo-url> --nodes 8 --no-submit

# Verbose output
python hpc_auto.py <repo-url> --nodes 16 --verbose
```

## Next Steps

1. **Run your test**: Start with a small test (4-8 nodes) to verify everything works
2. **Check output**: Monitor `output/<test_name>/` for job scripts and results
3. **Review results**: Check `efficiency_report.txt` for scaling analysis
4. **Scale up**: Once verified, run larger tests (16-64 nodes)
5. **Optimize**: Use custom launchers for better performance

## See Also

- `SYSTEM_CONFIG_GUIDE.md` - Detailed system configuration documentation
- `leonardo_system.py` - System configuration template
- `README.md` - Full framework documentation
