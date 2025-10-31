"""
Type definitions, enums, and constants for HPC-ScaleTest.
"""

from enum import Enum
from typing import Tuple


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


# Constants
DEFAULT_PROCS_PER_NODE = 128
DEFAULT_TIME_LIMIT = "01:00:00"
DEFAULT_OUTPUT_DIR = "output"
LARGE_JOB_THRESHOLD = 64  # Nodes threshold for special handling
