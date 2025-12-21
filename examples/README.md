# HPC-ScaleTest Examples

This directory contains example configurations and system definitions for HPC-ScaleTest.

## Files Overview

### Configuration Templates

- **`run.template.yaml`** - Comprehensive template with all available options
  - Use this as a starting point for your own configurations
  - Contains detailed comments explaining each setting
  - Copy and customize for your specific application

### Example Scaling Configurations

#### 2D Weak Scaling
- **File**: `example_2d_weak_scaling.yaml`
- **Use Case**: Applications where problem scales in 2 dimensions (e.g., thin slabs)
- **Pattern**: Scales X and Y dimensions while keeping Z constant
- **Progression**: 1 → 2 → 4 → 8 → 16 nodes

```bash
# Run 2D weak scaling test
python hpc_auto.py --config examples/example_2d_weak_scaling.yaml
```

#### 3D Weak Scaling
- **File**: `example_3d_weak_scaling.yaml`
- **Use Case**: Full 3D applications with uniform scaling
- **Pattern**: Scales X, Y, and Z dimensions proportionally
- **Progression**: 1 → 2 → 4 → 8 → 16 → 32 → 64 nodes

```bash
# Run 3D weak scaling test
python hpc_auto.py --config examples/example_3d_weak_scaling.yaml
```

#### Strong Scaling
- **File**: `example_strong_scaling.yaml`
- **Use Case**: Fixed problem size performance analysis
- **Pattern**: Constant problem size, increasing parallelism
- **Progression**: 1 → 2 → 4 → 8 → 16 → 32 nodes

```bash
# Run strong scaling test
python hpc_auto.py --config examples/example_strong_scaling.yaml
```

### System Configurations

- **`leonardo_system.py`** - Configuration for Leonardo supercomputer (CINECA)
  - Includes partition definitions (booster, dcgp)
  - Custom launcher configurations
  - Partition-specific resource settings

## Customizing Examples

### 1. Copy Template

```bash
cp examples/run.template.yaml my_application.yaml
```

### 2. Edit Configuration

Update these key sections for your application:

```yaml
# Your application repository
repository: https://github.com/username/my-hpc-app.git

# Problem size parameters
initial_procs: [4, 4, 2]        # 32 total processes
initial_domain: [100.0, 100.0, 50.0]
initial_cells: [256, 256, 128]

# Variable mapping (match your input file)
variable_map:
  length:
    x: "domain_size_x"  # Your variable name
    y: "domain_size_y"
    z: "domain_size_z"
```

### 3. Run Your Test

```bash
python hpc_auto.py --config my_application.yaml
```

## Common Customization Patterns

### Different Scaling Factors

Instead of doubling (factor 2.0), try other progressions:

```yaml
# Triple at each step: 1 → 3 → 9 → 27
scaling_factor: 3.0

# Or use specific node counts (overrides factor)
node_sequence: [1, 2, 4, 6, 8, 12, 16]
```

### GPU Applications

```yaml
hardware: gpu
gpus_per_node: 4
gpu_type: "A100"

# GPU-specific modules
modules:
  - cuda/11.8
  - gcc/11.2
  - openmpi/4.1
```

### Custom Build Flags

```yaml
build_system: cmake
build_flags:
  CMAKE_BUILD_TYPE: "Release"
  CMAKE_CXX_COMPILER: "g++"
  CMAKE_CXX_FLAGS: "-O3 -march=native -mtune=native"
  ENABLE_GPU: "ON"
```

### Multi-Tier QoS

```yaml
qos_mapping:
  debug:
    max_nodes: 2
    qos: "debug"
    time_limit: "00:30:00"
  normal:
    min_nodes: 3
    max_nodes: 16
    qos: "normal"
    time_limit: "02:00:00"
  large:
    min_nodes: 17
    max_nodes: 256
    qos: "large_scale"
    time_limit: "06:00:00"
```

## System Configuration Files

System configuration files (like `leonardo_system.py`) define:

1. **Available Partitions**: Hardware specifications for different queues
2. **Custom Launchers**: Optimized MPI launch configurations
3. **Module Environments**: System-specific software stacks
4. **Resource Policies**: QoS, time limits, and allocation rules

### Using System Configurations

```bash
# Command-line reference
python hpc_auto.py /path/to/code \
    --system-config examples/leonardo_system.py \
    --partition-name booster \
    --nodes 16

# Or in YAML
system_config: examples/leonardo_system.py
partition_name: booster
```

### Creating Your Own System Configuration

```python
# my_system.py
from utils.system_loader import SystemConfig, PartitionConfig

system = SystemConfig(
    name="MyCluster",
    scheduler="slurm",
    module_system="lmod",
)

# Define partition
system.add_partition(
    PartitionConfig(
        name="compute",
        nodes_available=256,
        procs_per_node=128,
        memory_per_node="256GB",
        has_gpu=False,
    )
)
```

## Best Practices

### 1. Start Small
- Test with `nodes: 4` first
- Verify input file parsing works
- Check job scripts before submission

### 2. Use Dry Run
```yaml
auto_submit: false  # Generate scripts only
```

### 3. Enable Verbose Logging
```yaml
verbose: true  # Detailed execution logs
```

### 4. Validate Before Running
```bash
# Check YAML syntax
python -c "import yaml; yaml.safe_load(open('my_config.yaml'))"

# Validate configuration
python hpc_auto.py --config my_config.yaml --dry-run
```

### 5. Keep Configurations in Version Control
```bash
git add my_application.yaml
git commit -m "Add scaling configuration for experiment X"
```

## Troubleshooting

### Configuration Not Found
```bash
# Use absolute path
python hpc_auto.py --config /full/path/to/config.yaml
```

### Variable Mapping Issues
- Check your input file variable names exactly
- Case-sensitive matching required
- Use `--verbose` to see parsing details

### Build Failures
```yaml
# Provide explicit build commands
build_system: make
build_flags:
  CC: "gcc"
  CXX: "g++"
```

## Additional Resources

- [Main README](../README.md) - Full framework documentation
- [YAML Configuration Guide](../YAML_CONFIG_GUIDE.md) - Complete YAML reference
- [Installation Guide](../INSTALL.md) - Setup instructions
- [Getting Started](../GETTING_STARTED.md) - Beginner's guide

## Questions?

- Review example configurations for similar use cases
- Check documentation in parent directory
- Open an issue on GitHub for help
