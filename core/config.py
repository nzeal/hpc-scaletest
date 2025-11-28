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
    qos_mapping: Optional[Dict[str, Dict]] = None  # QoS selection based on node count
    
    # Optional scheduler-specific settings
    partition: Optional[str] = None
    account: Optional[str] = None
    constraint: Optional[str] = None
    reservation: Optional[str] = None
    mail_user: Optional[str] = None
    
    def get_qos_for_nodes(self, num_nodes: int) -> Optional[str]:
        """
        Select appropriate QoS based on node count and QoS mapping.
        
        Args:
            num_nodes: Number of nodes for the job
            
        Returns:
            QoS string to use, or None if no mapping applies
            
        Example mapping:
            {
                "small": {"max_nodes": 16, "qos": "normal"},
                "large": {"min_nodes": 17, "qos": "dcgp_qos_bprod"}
            }
        """
        if not self.qos_mapping:
            # No mapping defined, return default QoS
            return self.qos
        
        # Find matching QoS tier based on node count
        for tier_name, tier_config in self.qos_mapping.items():
            max_nodes_allowed = tier_config.get('max_nodes', float('inf'))
            min_nodes_allowed = tier_config.get('min_nodes', 0)
            
            if min_nodes_allowed <= num_nodes <= max_nodes_allowed:
                selected_qos = tier_config.get('qos')
                if selected_qos:
                    return selected_qos
        
        # No matching tier found, return default
        return self.qos


@dataclass
class ScalingConfig:
    """Configuration for scaling tests."""
    scaling_type: ScalingType = ScalingType.STRONG
    max_nodes: int = 1
    
    # Initial decomposition and problem size
    initial_procs: ProcsDecomposition = (1, 1, 1)
    initial_domain: Optional[DomainSize] = None
    initial_cells: Optional[CellCount] = None
    particles_per_cell: Optional[ProcsDecomposition] = None  # (npcelx, npcely, npcelz)
    
    # Scaling factor (if defined, enables full weak scaling mode)
    scaling_factor: Optional[float] = None  # e.g., 2 (doubles per step)
    
    # Scaling dimensions: 2 for 2D (X→Y→X→Y), 3 for 3D (X→Y→Z→X→Y→Z)
    scaling_dimensions: int = 2  # Default to 2D scaling
    
    # Scaling factors (if defined, overrides node_sequence)
    weak_scaling_factors: Optional[List[float]] = None  # e.g., [1, 2, 4, 8]
    strong_scaling_factors: Optional[List[float]] = None  # e.g., [1, 2, 4, 8]
    
    # Node progression (default: powers of 2) - used only if scaling_factors not defined
    node_sequence: Optional[List[int]] = None
    
    # GENERIC: Variable mapping for input file parsing
    variable_map: Optional[Dict[str, Dict[str, str]]] = None  # e.g., {"length": {"x": "Lx", "y": "Ly", "z": "Lz"}}
    
    def get_scaling_factors(self) -> List[float]:
        """Get the scaling factors based on scaling type."""
        if self.scaling_type == ScalingType.WEAK:
            if self.weak_scaling_factors:
                return self.weak_scaling_factors
        elif self.scaling_type == ScalingType.STRONG:
            if self.strong_scaling_factors:
                return self.strong_scaling_factors
        
        # No scaling factors defined - return [1] for baseline mode
        return [1]
    
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
