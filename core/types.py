"""
Type definitions, enums, and constants for HPC-ScaleTest.
"""

import os
from enum import Enum
from typing import Tuple, Optional


class ScalingType(Enum):
    """Type of scaling test to perform."""
    STRONG = "strong"
    WEAK = "weak"


class SchedulerBackend(Enum):
    """Available job scheduler backends."""
    LOCAL = "local"
    SLURM = "slurm"
    PBS = "pbs"


class LauncherBackend(Enum):
    """Available MPI launcher backends."""
    LOCAL = "local"  # For login nodes / non-MPI
    SRUN = "srun"
    MPIRUN = "mpirun"
    MPIEXEC = "mpiexec"
    SIMPLE = "simple"


class ModuleBackend(Enum):
    """Available module system backends."""
    NOMOD = "nomod"
    TMOD = "tmod"
    TMOD4 = "tmod4"
    LMOD = "lmod"


class BuildBackend(Enum):
    """Available build system backends."""
    MAKE = "make"
    CMAKE = "cmake"
    AUTOTOOLS = "autotools"
    EASYBUILD = "easybuild"
    SPACK = "spack"


class JobStatus(Enum):
    """Job execution status."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"
    TIMEOUT = "timeout"
    UNKNOWN = "unknown"


# Type aliases for better code clarity
ProcsDecomposition = Tuple[int, int, int]  # (px, py, pz)
DomainSize = Tuple[float, float, float]    # (dx, dy, dz)
CellCount = Tuple[int, int, int]           # (nx, ny, nz)


def _detect_cpu_count() -> int:
    """Auto-detect number of CPUs on current system."""
    try:
        # Try to get physical cores from /proc/cpuinfo
        with open('/proc/cpuinfo', 'r') as f:
            content = f.read()
            # Count physical cores (not hyperthreads)
            cores = content.count('processor')
            if cores > 0:
                return cores
    except (FileNotFoundError, PermissionError):
        pass
    
    # Fallback to os.cpu_count()
    count = os.cpu_count()
    return count if count else 1


def get_default_procs_per_node() -> int:
    """
    Get default processes per node, auto-detecting if possible.
    
    Priority:
    1. HPC_SCALETEST_PROCS_PER_NODE environment variable
    2. Auto-detected CPU count
    3. Fallback to 1
    """
    env_value = os.environ.get('HPC_SCALETEST_PROCS_PER_NODE')
    if env_value:
        try:
            return int(env_value)
        except ValueError:
            pass
    
    return _detect_cpu_count()


# Constants - use functions for dynamic detection where possible
DEFAULT_PROCS_PER_NODE = get_default_procs_per_node()
DEFAULT_TIME_LIMIT = os.environ.get('HPC_SCALETEST_TIME_LIMIT', "01:00:00")
DEFAULT_OUTPUT_DIR = os.environ.get('HPC_SCALETEST_OUTPUT_DIR', "output")
LARGE_JOB_THRESHOLD = int(os.environ.get('HPC_SCALETEST_LARGE_JOB_THRESHOLD', "64"))
