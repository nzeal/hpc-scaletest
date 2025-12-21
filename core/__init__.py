# ====================
# core/__init__.py
# ====================
"""
Core abstractions and types for HPC-ScaleTest.

Key modules:
- unified_execution: SINGLE SOURCE OF TRUTH for execution logic (NEW)
- hardware: Unified hardware topology detection
- topology: Centralized hardware topology detection
- mpi_command: MPI command generation
- job_script: SLURM job script generation
- slurm_script: Clean SLURM script generator
- result_writer: Incremental result writer
"""

from .types import (
    ScalingType,
    SchedulerBackend,
    LauncherBackend,
    ModuleBackend,
    BuildBackend,
    JobStatus
)
from .test_definition import Test
from .factory import BackendFactory
from .registry import (
    register_launcher,
    get_launcher,
    list_launchers,
    has_launcher,
    JobLauncher
)

# Import new UNIFIED EXECUTION module (single source of truth)
try:
    from .unified_execution import (
        # Data classes
        HardwareTopology as UnifiedTopology,
        MPIConfiguration as UnifiedMPIConfig,
        SlurmConfiguration as UnifiedSlurmConfig,
        RunResult as UnifiedRunResult,
        # Enums
        GPUVendor as UnifiedGPUVendor,
        ScalingType as UnifiedScalingType,
        # Functions
        detect_topology as unified_detect_topology,
        generate_gpu_binding_script as unified_gpu_bind_script,
        generate_job_script as unified_job_script,
        generate_node_sequence,
        get_detector,
        # Classes
        TopologyDetector as UnifiedTopologyDetector,
        IncrementalResultWriter as UnifiedResultWriter,
        # Constants
        GPU_BIND_SCRIPT,
    )
    HAS_UNIFIED_EXECUTION = True
except ImportError:
    HAS_UNIFIED_EXECUTION = False

# Import new unified hardware module
try:
    from .hardware import (
        HardwareTopology,
        MPIConfiguration,
        TopologyDetector as HardwareTopologyDetector,
        detect_topology as detect_hardware,
        get_mpi_config,
        generate_gpu_bind_script,
    )
    HAS_HARDWARE = True
except ImportError:
    HAS_HARDWARE = False

# Import new SLURM script module
try:
    from .slurm_script import (
        SlurmJobConfig,
        generate_slurm_job_script,
        create_job_config_from_topology,
    )
    HAS_SLURM_SCRIPT = True
except ImportError:
    HAS_SLURM_SCRIPT = False

# Import new result writer module
try:
    from .result_writer import (
        RunResult,
        ScalingSummary,
        IncrementalResultWriter,
        create_result_writer,
    )
    HAS_RESULT_WRITER = True
except ImportError:
    HAS_RESULT_WRITER = False

# Import topology module (central hardware detection)
try:
    from .topology import (
        NodeTopology,
        MPIMapping,
        GPUVendor,
        DetectionContext,
        TopologyDetector,
        get_topology_detector,
        detect_topology,
        compute_mpi_mapping,
    )
    HAS_TOPOLOGY = True
except ImportError:
    HAS_TOPOLOGY = False

# Import MPI command module
try:
    from .mpi_command import (
        MPIImplementation,
        MPIInfo,
        MPIDetector,
        MPICommandGenerator,
        generate_gpu_binding_script,
        get_mpi_command,
    )
    HAS_MPI_COMMAND = True
except ImportError:
    HAS_MPI_COMMAND = False

# Import job script generator
try:
    from .job_script import (
        JobConfiguration,
        JobScriptGenerator,
        generate_job_script,
    )
    HAS_JOB_SCRIPT = True
except ImportError:
    HAS_JOB_SCRIPT = False

__all__ = [
    # Types
    'ScalingType',
    'SchedulerBackend',
    'LauncherBackend',
    'ModuleBackend',
    'BuildBackend',
    'JobStatus',
    # Core classes
    'Test',
    'BackendFactory',
    # Registry
    'register_launcher',
    'get_launcher',
    'list_launchers',
    'has_launcher',
    'JobLauncher',
]

# Add unified execution exports if available (RECOMMENDED)
if HAS_UNIFIED_EXECUTION:
    __all__.extend([
        'UnifiedTopology',
        'UnifiedMPIConfig',
        'UnifiedSlurmConfig',
        'UnifiedRunResult',
        'UnifiedGPUVendor',
        'UnifiedScalingType',
        'unified_detect_topology',
        'unified_gpu_bind_script',
        'unified_job_script',
        'generate_node_sequence',
        'get_detector',
        'UnifiedTopologyDetector',
        'UnifiedResultWriter',
        'GPU_BIND_SCRIPT',
    ])

# Add hardware exports if available
if HAS_HARDWARE:
    __all__.extend([
        'HardwareTopology',
        'MPIConfiguration',
        'HardwareTopologyDetector',
        'detect_hardware',
        'get_mpi_config',
        'generate_gpu_bind_script',
    ])

# Add SLURM script exports if available
if HAS_SLURM_SCRIPT:
    __all__.extend([
        'SlurmJobConfig',
        'generate_slurm_job_script',
        'create_job_config_from_topology',
    ])

# Add result writer exports if available
if HAS_RESULT_WRITER:
    __all__.extend([
        'RunResult',
        'ScalingSummary',
        'IncrementalResultWriter',
        'create_result_writer',
    ])

# Add topology exports if available
if HAS_TOPOLOGY:
    __all__.extend([
        'NodeTopology',
        'MPIMapping',
        'GPUVendor',
        'DetectionContext',
        'TopologyDetector',
        'get_topology_detector',
        'detect_topology',
        'compute_mpi_mapping',
    ])

# Add MPI command exports if available
if HAS_MPI_COMMAND:
    __all__.extend([
        'MPIImplementation',
        'MPIInfo',
        'MPIDetector',
        'MPICommandGenerator',
        'generate_gpu_binding_script',
        'get_mpi_command',
    ])

# Add job script exports if available
if HAS_JOB_SCRIPT:
    __all__.extend([
        'JobConfiguration',
        'JobScriptGenerator',
        'generate_job_script',
    ])
