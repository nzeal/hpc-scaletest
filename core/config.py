"""
Configuration dataclasses for test definition.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, List

from .types import (
    ScalingType, SchedulerBackend, LauncherBackend, 
    ModuleBackend, BuildBackend,
    ProcsDecomposition, DomainSize, CellCount,
    DEFAULT_PROCS_PER_NODE, DEFAULT_TIME_LIMIT
)


@dataclass
class BackendConfig:
    """Configuration for backend systems."""
    scheduler: SchedulerBackend = SchedulerBackend.LOCAL
    launcher: LauncherBackend = LauncherBackend.SRUN
    module_system: ModuleBackend = ModuleBackend.NOMOD
    build_system: Optional[BuildBackend] = None
    
    # Backend-specific options
    scheduler_options: Dict = field(default_factory=dict)
    launcher_options: Dict = field(default_factory=dict)
    module_options: Dict = field(default_factory=dict)
    build_options: Dict = field(default_factory=dict)


@dataclass
class ResourceConfig:
    """Configuration for computational resources."""
    max_nodes: int = 1
    procs_per_node: int = DEFAULT_PROCS_PER_NODE
    gpus_per_node: int = 0
    memory_per_node: Optional[str] = None
    time_limit: str = DEFAULT_TIME_LIMIT
    exclusive: bool = False
    qos: Optional[str] = None
    
    # Optional scheduler-specific settings
    partition: Optional[str] = None
    account: Optional[str] = None
    constraint: Optional[str] = None
    reservation: Optional[str] = None
    mail_user: Optional[str] = None


@dataclass
class ScalingConfig:
    """Configuration for scaling tests."""
    scaling_type: ScalingType = ScalingType.STRONG
    max_nodes: int = 1
    
    # Initial decomposition and problem size
    initial_procs: ProcsDecomposition = (1, 1, 1)
    initial_domain: Optional[DomainSize] = None
    initial_cells: Optional[CellCount] = None
    
    # Node progression (default: powers of 2)
    node_sequence: Optional[List[int]] = None
    
    def get_node_sequence(self) -> List[int]:
        """Get the sequence of node counts for scaling."""
        if self.node_sequence:
            return self.node_sequence
        
        # Generate powers of 2 up to max_nodes
        sequence = []
        n = 1
        while n <= self.max_nodes:
            sequence.append(n)
            n *= 2
        
        # Ensure max_nodes is included if not a power of 2
        if sequence[-1] != self.max_nodes:
            sequence.append(self.max_nodes)
        
        return sequence


@dataclass
class EnvironmentConfig:
    """Configuration for environment setup."""
    modules: List[str] = field(default_factory=list)
    env_vars: Dict[str, str] = field(default_factory=dict)
    pre_commands: List[str] = field(default_factory=list)
    post_commands: List[str] = field(default_factory=list)


@dataclass
class BuildConfig:
    """Configuration for building the application."""
    source_dir: Optional[Path] = None
    build_dir: Optional[Path] = None
    install_dir: Optional[Path] = None
    build_flags: Dict[str, str] = field(default_factory=dict)
    parallel_jobs: int = 4
    executable_name: str = "main_exe"  # e.g., "iPIC3D" or "ECsim"


@dataclass
class JobConfig:
    """Configuration for a single job instance."""
    job_id: str
    num_nodes: int
    num_procs: int
    procs_decomposition: ProcsDecomposition
    
    # Problem-specific parameters
    domain_size: Optional[DomainSize] = None
    cell_count: Optional[CellCount] = None
    
    # Runtime parameters
    working_dir: Optional[Path] = None
    output_file: Optional[Path] = None
    
    def total_procs(self) -> int:
        """Calculate total number of processes."""
        px, py, pz = self.procs_decomposition
        return px * py * pz
