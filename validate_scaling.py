"""Quick validation test for new scaling engine."""

import sys
sys.path.insert(0, '.')

from engine.scaling import ScalingEngine
from core.config import ScalingConfig, ResourceConfig
from core.types import ScalingType

# Test 2D scaling
config = ScalingConfig(
    scaling_type=ScalingType.WEAK,
    max_nodes=8,
    scaling_factor=2.0,
    scaling_dimensions=2,
    initial_procs=(14, 8, 1),
    initial_domain=(84.0, 48.0, 1.0),
    initial_cells=(840, 480, 1),
    particles_per_cell=(20, 20, 1)
)

resource_config = ResourceConfig(procs_per_node=112)
engine = ScalingEngine(config, resource_config)

jobs = engine.generate_job_configs()

print("\n" + "="*80)
print("VALIDATION TEST: 2D Weak Scaling")
print("="*80)

expected = {
    1: {"domain": (84.0, 48.0, 1.0), "cells": (840, 480, 1), "procs": (14, 8, 1)},
    2: {"domain": (168.0, 48.0, 1.0), "cells": (1680, 480, 1), "procs": (28, 8, 1)},
    4: {"domain": (168.0, 96.0, 1.0), "cells": (1680, 960, 1), "procs": (28, 16, 1)},
    8: {"domain": (336.0, 96.0, 1.0), "cells": (3360, 960, 1), "procs": (56, 16, 1)}
}

all_pass = True
for job in jobs:
    exp = expected[job.num_nodes]
    match = (
        job.domain_size == exp["domain"] and
        job.cell_count == exp["cells"] and
        job.procs_decomposition == exp["procs"]
    )
    
    status = "\u2713 PASS" if match else "\u2717 FAIL"
    print(f"\nNode {job.num_nodes}: {status}")
    print(f"  Domain:     {job.domain_size} (expected {exp['domain']})")
    print(f"  Cells:      {job.cell_count} (expected {exp['cells']})")
    print(f"  MPI decomp: {job.procs_decomposition} (expected {exp['procs']})")
    
    if not match:
        all_pass = False

print("\n" + "="*80)
if all_pass:
    print("\u2705 ALL TESTS PASSED!")
else:
    print("\u274c SOME TESTS FAILED")
print("="*80)
