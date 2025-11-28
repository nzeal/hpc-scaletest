# Scaling Engine Integration Guide

## Overview

The new scaling engine provides deterministic, rule-based weak scaling for HPC applications with full support for 1D, 2D, and 3D modes.

## Files Created

1. **`engine/scaling.py`** - Production scaling engine
2. **`engine/scaling_debug.py`** - Debug-enabled version with verbose logging
3. **`tests/test_scaling_engine.py`** - Comprehensive pytest test suite

## Key Features

✅ **Deterministic scaling** - Same inputs always produce same outputs  
✅ **Dimension-aware** - Supports 1D (X only), 2D (X→Y), 3D (X→Y→Z) patterns  
✅ **run.yaml driven** - All parameters from configuration, no hardcoded physics  
✅ **Node 1 baseline override** - Exact run.yaml values, ignores input file  
✅ **MPI rank correction** - Automatic adjustment to match `nodes × procs_per_node`  
✅ **Constant particles/cell** - Maintains particle density across all nodes  

## Scaling Rules

### 1D Mode (`scaling_dimensions: 1`)
- **Pattern**: X only
- **Constant**: Y, Z never change
- **Example**:
  ```
  Node 1: (Lx, Ly, Lz) = (84, 48, 1)
  Node 2: (168, 48, 1)  ← X×2
  Node 4: (336, 48, 1)  ← X×4
  ```

### 2D Mode (`scaling_dimensions: 2`)
- **Pattern**: X → Y → X → Y → ...
- **Constant**: Z never changes
- **Example**:
  ```
  Node 1: (Lx, Ly, Lz) = (84, 48, 1)
  Node 2: (168, 48, 1)  ← X×2
  Node 4: (168, 96, 1)  ← Y×2
  Node 8: (336, 96, 1)  ← X×4
  ```

### 3D Mode (`scaling_dimensions: 3`)
- **Pattern**: X → Y → Z → X → Y → Z → ...
- **Example**:
  ```
  Node 1: (Lx, Ly, Lz) = (10, 10, 10)
  Node 2: (20, 10, 10)  ← X×2
  Node 4: (20, 20, 10)  ← Y×2
  Node 8: (20, 20, 20)  ← Z×2
  ```

## Usage

### Standard Mode

```python
from engine.scaling import ScalingEngine
from core.config import ScalingConfig, ResourceConfig
from core.types import ScalingType

# Configuration from run.yaml
config = ScalingConfig(
    scaling_type=ScalingType.WEAK,
    max_nodes=8,
    scaling_factor=2.0,
    scaling_dimensions=2,  # 2D scaling
    initial_procs=(14, 8, 1),
    initial_domain=(84.0, 48.0, 1.0),
    initial_cells=(840, 480, 1),
    particles_per_cell=(20, 20, 1)
)

resource_config = ResourceConfig(procs_per_node=112)

# Create engine
engine = ScalingEngine(config, resource_config)

# Generate job configurations
jobs = engine.generate_job_configs()

# Use jobs for submission
for job in jobs:
    print(f"{job.job_id}: {job.num_nodes} nodes, decomp={job.procs_decomposition}")
```

### Debug Mode

```python
from engine.scaling_debug import DebugScalingEngine

# Create debug engine with verbose output
engine = DebugScalingEngine(
    config,
    resource_config,
    debug_level='VERBOSE'  # Options: 'minimal', 'standard', 'verbose', 'trace'
)

jobs = engine.generate_job_configs()
```

Debug levels:
- `minimal` - Only final results
- `standard` - Key steps and corrections
- `verbose` - Full step-by-step details
- `trace` - Every calculation logged

## Run.yaml Configuration

Required fields for weak scaling:

```yaml
# Scaling configuration
scaling: weak
nodes: 8
scaling_factor: 2.0
scaling_dimensions: 2  # 1, 2, or 3

# Baseline parameters (Node 1)
initial_procs: [14, 8, 1]       # px, py, pz
initial_domain: [84.0, 48.0, 1.0]  # Lx, Ly, Lz
initial_cells: [840, 480, 1]    # nx, ny, nz
particles_per_cell: {x: 20, y: 20, z: 1}

# Hardware
procs_per_node: 112

# Variable mapping (tells parser which input file variables to modify)
variable_map:
  length: {x: "Lx", y: "Ly", z: "Lz"}
  cells: {x: "nxc", y: "nyc", z: "nzc"}
  mpi_decomposition: {x: "XLEN", y: "YLEN", z: "ZLEN"}
  particles: {x: "npcelx", y: "npcely", z: "npcelz"}
```

## Testing

Run the comprehensive test suite:

```bash
# Run all tests
pytest tests/test_scaling_engine.py -v

# Run specific test class
pytest tests/test_scaling_engine.py::TestWeakScaling2D -v

# Run with coverage
pytest tests/test_scaling_engine.py --cov=engine.scaling --cov-report=html
```

Test coverage includes:
- ✅ 1D scaling pattern validation
- ✅ 2D alternating X→Y pattern
- ✅ 3D cyclic X→Y→Z pattern
- ✅ MPI rank correction
- ✅ Node 1 baseline override
- ✅ Constant particles per cell
- ✅ Configuration validation
- ✅ Error handling

## Migration from Old Engine

The new engine is a drop-in replacement. Simply use:

```python
from engine.scaling import ScalingEngine  # New deterministic engine
```

instead of the old implementation.

### Key Differences

| Feature | Old Engine | New Engine |
|---------|-----------|------------|
| Scaling logic | Incremental (previous step) | Cumulative (from baseline) |
| Node 1 | Sometimes modified | Always exact run.yaml |
| MPI correction | Multiple dimensions | Active dimension only |
| Debugging | Limited logging | Full debug mode available |
| Testing | Minimal | Comprehensive pytest suite |

## Validation

To verify correctness:

1. **Run debug mode**:
   ```python
   engine = DebugScalingEngine(config, resource_config, debug_level='VERBOSE')
   jobs = engine.generate_job_configs()
   ```

2. **Check logs** - Verify:
   - Active dimension matches expected pattern
   - Multipliers are correct (`factor^count`)
   - MPI decomposition matches `nodes × procs_per_node`
   - Z is constant in 2D mode

3. **Run tests**:
   ```bash
   pytest tests/test_scaling_engine.py -v
   ```

## Troubleshooting

### Issue: MPI ranks don't match

**Symptoms**: Error `"MPI correction failed"`

**Solution**: Check that `procs_per_node` in run.yaml matches your hardware and that the scaling factor produces valid decompositions.

### Issue: Z dimension changes in 2D mode

**Symptoms**: `pz`, `nz`, or `Lz` differ from baseline

**Solution**: Ensure `scaling_dimensions: 2` is set in run.yaml. The engine enforces Z constancy in 2D mode.

### Issue: Node 1 doesn't match run.yaml

**Symptoms**: Baseline parameters differ from configuration

**Solution**: The new engine always uses exact run.yaml values for Node 1. Verify your run.yaml has the correct baseline.

## Example Output

With `debug_level='VERBOSE'`, you'll see:

```
╔══════════════════════════════════════════════════════════════════════════════╗
║ WEAK SCALING - 2D MODE (factor=2.0)                                          ║
║ Baseline (run.yaml): L=(84.0, 48.0, 1.0)                                     ║
║                      n=(840, 480, 1)                                          ║
║                      p=(14, 8, 1)                                             ║
╚══════════════════════════════════════════════════════════════════════════════╝

================================================================================
NODE 1: BASELINE (run.yaml override)
================================================================================
  Domain:     (84.0, 48.0, 1.0)
  Cells:      (840, 480, 1)
  MPI decomp: (14, 8, 1) = 112 procs
  Particles:  (20, 20, 1)
  ⚠ Input file values IGNORED - using run.yaml

################################################################################
# STEP 1 → NODE 2
################################################################################

┌─ DIMENSION SELECTION ─────────────────────────────────────
│ Scaling mode: 2D
│ Step index:   1
│ Active dim:   X
└────────────────────────────────────────────────────────────

┌─ CUMULATIVE SCALE COUNTS ─────────────────────────────────
│ X has scaled 1 times → multiplier = 2.0^1 = 2.000
│ Y has scaled 0 times → multiplier = 2.0^0 = 1.000
│ Z has scaled 0 times → multiplier = 2.0^0 = 1.000
└────────────────────────────────────────────────────────────

┌─ SCALED VALUES (before MPI correction) ───────────────────
│ Domain: Lx=168.00 (84.00×2.00)
│         Ly=48.00 (48.00×1.00)
│         Lz=1.00 (1.00×1.00)
│
│ Cells:  nx=1680 (840×2.00)
│         ny=480 (480×1.00)
│         nz=1 (1×1.00)
│
│ Procs:  px=28 (14×2.00)
│         py=8 (8×1.00)
│         pz=1 (1×1.00)
└────────────────────────────────────────────────────────────

┌─ MPI RANK VALIDATION ─────────────────────────────────────
│ Required: 2 nodes × 112 procs/node = 224
│ Computed: 28×8×1 = 224
│ ✓ No correction needed
└────────────────────────────────────────────────────────────

┌─ FINAL CONFIGURATION ─────────────────────────────────────
│ Node:       2
│ Domain:     (168.00, 48.00, 1.00)
│ Cells:      (1680, 480, 1)
│ MPI decomp: (28, 8, 1) = 224
│ Particles:  (20, 20, 1) [CONSTANT]
└────────────────────────────────────────────────────────────
```

## Support

For issues or questions:
1. Check the test suite for examples
2. Enable `debug_level='VERBOSE'` to see detailed logs
3. Verify run.yaml configuration
4. Check that `variable_map` matches your input file format
