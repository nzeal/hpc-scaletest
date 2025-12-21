# ====================
# engine/__init__.py
# ====================
"""
Test execution engine for HPC-ScaleTest.

Modules:
- scaling: Scaling configuration generation (weak/strong)
- runner: Test execution and job submission
- cpu_execution: CPU-only job execution
- gpu_execution: GPU-accelerated job execution
"""

from .scaling import ScalingEngine
from .runner import TestRunner

# Import new execution modules
try:
    from .cpu_execution import CPUExecutionEngine, CPUJobConfig
    from .gpu_execution import GPUExecutionEngine, GPUJobConfig
    HAS_EXECUTION_MODULES = True
except ImportError:
    HAS_EXECUTION_MODULES = False

__all__ = [
    'ScalingEngine',
    'TestRunner',
    'CPUExecutionEngine',
    'CPUJobConfig',
    'GPUExecutionEngine', 
    'GPUJobConfig',
]
