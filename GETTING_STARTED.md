# Getting Started with HPC Auto

Welcome to HPC Auto - the automated end-to-end framework for HPC scaling tests!

---

## Installation

No additional installation required! The framework uses Python standard library and existing dependencies.

```bash
cd hpc-scaletest2
# You're ready to go!
```

---

## 5-Minute Quick Start

### 1. Run the Demo (No Actual Tests)

```bash
python demo_auto_workflow.py
```

This shows you all capabilities without running real tests.

---

### 2. Your First Real Test

**Option A: With Local Code**
```bash
python hpc_auto.py /path/to/your/hpc/code \
    --scaling strong \
    --nodes 2 \
    --no-submit \
    --verbose
```

**Option B: From Git Repository**
```bash
python hpc_auto.py https://github.com/user/hpc-app.git \
    --scaling strong \
    --nodes 2 \
    --no-submit \
    --verbose
```

The `--no-submit` flag generates job scripts without submitting them, so you can review first.

---

### 3. Review Generated Files

```bash
# Check the output
ls output/*/

# Look at a job script
cat output/*/nodes_1/job.sh

# Review detected configuration
cat output/*/nodes_1/metadata.json
```

---

### 4. Submit Jobs Manually (if --no-submit was used)

```bash
# Submit each job
sbatch output/*/nodes_1/job.sh
sbatch output/*/nodes_2/job.sh

# Or use the batch submitter
python utils/job_submitter.py output/your_test_dir/
```

---

### 5. View Reports (After Jobs Complete)

```bash
# Wait for jobs to complete, then:
cat output/*/StrongScalingReport.txt
```

---

## Common Use Cases

### Use Case 1: Quick Strong Scaling Test

**Scenario**: You have a code and want to test strong scaling up to 8 nodes.

```bash
python hpc_auto.py /path/to/code --scaling strong --nodes 8
```

**What happens:**
- ✅ Code analyzed for dependencies
- ✅ Automatically compiled
- ✅ Tests generated for 1, 2, 4, 8 nodes
- ✅ Jobs submitted automatically
- ✅ Efficiency report generated

---

### Use Case 2: Weak Scaling with Custom Problem Size

**Scenario**: Test how your code scales when problem size grows proportionally.

```bash
python hpc_auto.py /path/to/code \
    --scaling weak \
    --nodes 16 \
    --initial-domain "10.0,10.0,10.0" \
    --initial-cells "256,256,256"
```

---

### Use Case 3: GPU-Enabled Application

**Scenario**: Test a CUDA application on GPU nodes.

```bash
python hpc_auto.py /path/to/cuda-app \
    --hardware gpu \
    --nodes 4 \
    --procs-per-node 4 \
    --gpus-per-node 1 \
    --partition gpu_prod
```

---

### Use Case 4: Custom Modules and Configuration

**Scenario**: Your system has specific module requirements.

```bash
python hpc_auto.py /path/to/code \
    --nodes 8 \
    --modules "gcc/11.2.0,openmpi/4.1.1,hdf5/1.12.0" \
    --build-system cmake \
    --partition production \
    --account my_project
```

---

### Use Case 5: Generate Scripts for Manual Review

**Scenario**: You want to review job scripts before submission.

```bash
python hpc_auto.py /path/to/code \
    --nodes 4 \
    --no-submit \
    --verbose

# Review scripts
cat output/*/nodes_*/job.sh

# Submit when satisfied
cd output/your_test_dir/
for dir in nodes_*/; do
    sbatch $dir/job.sh
done
```

---

## Python API Usage

### Simple API

```python
from engine.orchestrator import create_simple_workflow
from utils.logging_config import setup_logging

# Enable logging
setup_logging()

# Create workflow
orchestrator = create_simple_workflow(
    source="/path/to/code",
    scaling_type="strong",
    max_nodes=8,
    hardware_type="cpu"
)

# Run it
success = orchestrator.run()

if success:
    print("✅ Test completed successfully!")
else:
    print("❌ Test failed")
```

### Advanced API

```python
from engine.orchestrator import HPCOrchestrator, OrchestratorConfig
from pathlib import Path

# Detailed configuration
config = OrchestratorConfig(
    # Required
    source="/path/to/code",
    
    # Scaling
    scaling_type="strong",
    max_nodes=16,
    initial_procs=(4, 4, 2),
    
    # Resources
    procs_per_node=128,
    time_limit="04:00:00",
    partition="production",
    account="my_project",
    
    # Build
    modules=["gcc/11.2.0", "openmpi/4.1.1"],
    build_system="cmake",
    
    # Output
    output_dir=Path("my_results"),
    test_name="my_benchmark"
)

# Create and run
orchestrator = HPCOrchestrator(config)
success = orchestrator.run()
```

---

## Understanding the Output

### Directory Structure

```
output/
└── my_app_strong_20251025_104500/
    ├── nodes_1/              # First test (1 node)
    │   ├── job.sh           # Submission script
    │   ├── input.dat        # Scaled input
    │   ├── metadata.json    # Job metadata
    │   └── out_*.out        # Output (after run)
    │
    ├── nodes_2/              # Second test (2 nodes)
    ├── nodes_4/              # Third test (4 nodes)
    ├── nodes_8/              # Fourth test (8 nodes)
    │
    ├── install/              # Compiled code
    │   └── my_app
    │
    ├── summary.json                # Test summary
    ├── StrongScalingReport.txt     # Efficiency report
    └── strong_scaling_report.json  # JSON data
```

### Example Report

```
==============================================================================
Strong Scaling Efficiency Report
==============================================================================
Test Name: my_app
==============================================================================

Nodes      Procs      Time(s)      Speedup      Efficiency
------------------------------------------------------------------------------
1          128        100.00       1.00         100.0%
2          256        52.00        1.92         96.2%
4          512        28.00        3.57         89.3%
8          1024       16.00        6.25         78.1%
==============================================================================

Summary Statistics:
  Average Efficiency: 90.9%
  Maximum Speedup: 6.25x
==============================================================================
```

**Interpretation:**
- **Speedup**: How much faster compared to 1 node (ideal: linear)
- **Efficiency**: Resource utilization percentage (ideal: 100%)
- **Good scaling**: Efficiency > 70% is typically considered good

---

## Tips for Success

### 1. **Improve Auto-Detection**

Add good documentation to your code:

```markdown
# README.md

## Dependencies
- MPI (OpenMPI 4.1+)
- HDF5 1.12+
- CUDA 11.0+ (optional)

## Building
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
```
```

### 2. **Start Small**

Always test with small node counts first:
```bash
python hpc_auto.py /path/to/code --nodes 2 --verbose
```

### 3. **Use --verbose or --debug**

See what's happening:
```bash
python hpc_auto.py /path/to/code --nodes 4 --verbose
```

### 4. **Test Locally First**

Use local scheduler for development:
```bash
python hpc_auto.py /path/to/code --scheduler local --nodes 1
```

---

## Troubleshooting

### Problem: "Git clone failed"

**Solution:**
```bash
# Check Git is installed
git --version

# Try cloning manually first
git clone https://your-repo-url.git
```

---

### Problem: "Build system not detected"

**Solution:**
```bash
# Manually specify
python hpc_auto.py /path/to/code --build-system cmake
```

---

### Problem: "Modules not found"

**Solution:**
```bash
# Check available modules
module avail

# Specify exact names
python hpc_auto.py /path/to/code \
    --modules "gcc/11.2.0,openmpi/4.1.1"
```

---

### Problem: "No timing data in report"

**Solution:** Ensure your application outputs timing in a recognizable format:
- "Time: 123.45"
- "Elapsed: 123.45 seconds"
- "Wall time: 123.45"

---

## Next Steps

1. **Read Detailed Docs**
   - [Automated Workflow Guide](docs/AUTOMATED_WORKFLOW.md) - Complete reference
   - [Quick Start Guide](docs/QUICKSTART.md) - More examples

2. **Explore Examples**
   - `examples/simple_workflow.py` - Simple API
   - `examples/advanced_workflow.py` - Advanced configuration

3. **Get Help**
   ```bash
   python hpc_auto.py --help
   ```

4. **Review Implementation**
   - [Implementation Summary](IMPLEMENTATION_SUMMARY.md) - Technical details

---

## Command Reference

### Most Common Commands

```bash
# Basic usage
python hpc_auto.py /path/to/code --nodes 8

# With options
python hpc_auto.py /path/to/code \
    --scaling strong \
    --nodes 16 \
    --partition production \
    --account my_project \
    --verbose

# From Git
python hpc_auto.py https://github.com/user/repo.git --nodes 8

# GPU test
python hpc_auto.py /path/to/code --hardware gpu --nodes 4

# Generate only (no submit)
python hpc_auto.py /path/to/code --nodes 4 --no-submit
```

### Get All Options

```bash
python hpc_auto.py --help
```

---

## Support

For issues or questions:
1. Check troubleshooting section above
2. Review documentation in `docs/`
3. Run with `--verbose` or `--debug` for detailed logs
4. Check generated job scripts in `output/`

---

## Summary

**What You Need:**
- Python 3.7+
- HPC application source code (local or Git)
- Access to HPC cluster (or local machine for testing)

**What You Get:**
- ✅ Automated dependency detection
- ✅ Automatic compilation
- ✅ Scaling test generation
- ✅ Job submission and monitoring
- ✅ Beautiful efficiency reports

**How to Start:**
```bash
python hpc_auto.py /path/to/code --nodes 8
```

**That's it!** 🚀

---

**Happy scaling!**
