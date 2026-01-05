"""
Base Job Generator Module

Provides abstract base class and common configuration for job generators.
All hardware-specific generators inherit from JobGeneratorBase.

Design Principles:
==================
- NO HARDCODED VALUES: All values from configuration
- SYSTEM AGNOSTIC: Works on any HPC system
- CONFIGURABLE: All parameters passed via JobGeneratorConfig

Author: HPC-ScaleTest Contributors
"""

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class JobGeneratorConfig:
    """
    Configuration for job script generation.
    
    IMPORTANT: All values MUST be provided by the user via YAML configuration.
    There are NO hardcoded defaults - this ensures system and application agnostic behavior.
    
    Attributes:
        # Hardware Configuration (REQUIRED - from YAML)
        cpus_per_node: CPU cores per node (from system detection or YAML)
        gpus_per_node: GPUs per node (0 for CPU jobs)
        
        # Resource Configuration (REQUIRED - from YAML)
        partition: SLURM partition name
        account: SLURM account/project
        qos: Quality of service
        time_limit: Job time limit (HH:MM:SS)
        exclusive: Request exclusive node access
        
        # Application Configuration (REQUIRED - from YAML)
        binary_dir: Directory containing the compiled binary
        executable: Name of the executable
        input_file: Input file name
        
        # MPI Configuration
        launcher: MPI launcher ('mpirun', 'srun')
        procs_per_node: MPI processes per node (calculated for GPU jobs)
        
        # Environment Configuration
        modules: List of modules to load
        env_vars: Environment variables to set
        
        # Scaling Configuration (REQUIRED - from YAML)
        max_nodes: Maximum nodes for scaling tests
        node_sequence: Explicit node sequence (if provided)
    """
    # Hardware - NO DEFAULTS, must be provided
    cpus_per_node: int
    gpus_per_node: int = 0
    
    # Resources - NO DEFAULTS for critical fields
    partition: str = ""
    account: str = ""
    qos: str = ""
    time_limit: str = "02:00:00"
    exclusive: bool = True
    
    # Application - NO DEFAULTS, must be provided
    binary_dir: str = ""
    executable: str = ""
    input_file: str = ""
    
    # MPI
    launcher: str = "mpirun"
    procs_per_node: int = 0  # Calculated based on hardware type
    
    # Environment
    modules: List[str] = field(default_factory=list)
    env_vars: Dict[str, str] = field(default_factory=dict)
    
    # Scaling - NO DEFAULTS, must be provided
    max_nodes: int = 1
    node_sequence: Optional[List[int]] = None
    
    def get_node_sequence(self) -> List[int]:
        """
        Get the sequence of node counts for scaling.
        
        If node_sequence is provided, use it directly.
        Otherwise, generate powers of 2 up to max_nodes.
        
        Returns:
            List of node counts
        """
        if self.node_sequence:
            return self.node_sequence
        
        # Generate powers of 2 up to max_nodes
        sequence = []
        n = 1
        while n <= self.max_nodes:
            sequence.append(n)
            n *= 2
        
        # Ensure max_nodes is included if not a power of 2
        if sequence and sequence[-1] != self.max_nodes:
            sequence.append(self.max_nodes)
        
        return sequence


class JobGeneratorBase(ABC):
    """
    Abstract base class for job script generators.
    
    Subclasses must implement:
    - generate_mpi_command(): Generate the MPI launch command
    - generate_script(): Generate the complete job script
    
    Common functionality provided:
    - SLURM directive generation
    - Environment setup
    - Working directory handling
    """
    
    def __init__(self, config: JobGeneratorConfig):
        """
        Initialize job generator with configuration.
        
        Args:
            config: Job generator configuration
        """
        self.config = config
    
    @abstractmethod
    def generate_mpi_command(
        self,
        num_nodes: int,
        total_ranks: int,
        ranks_per_node: int,
        cores_per_rank: int
    ) -> str:
        """
        Generate the MPI launch command.
        
        Args:
            num_nodes: Number of nodes
            total_ranks: Total MPI ranks
            ranks_per_node: MPI ranks per node
            cores_per_rank: CPU cores per rank
        
        Returns:
            MPI launch command as string
        """
        pass
    
    @abstractmethod
    def generate_script(
        self,
        num_nodes: int,
        job_name: str,
        working_dir: str,
        **kwargs
    ) -> str:
        """
        Generate complete job script.
        
        Args:
            num_nodes: Number of nodes
            job_name: SLURM job name
            working_dir: Working directory
            **kwargs: Additional parameters
        
        Returns:
            Complete job script as string
        """
        pass
    
    def generate_slurm_directives(
        self,
        num_nodes: int,
        job_name: str,
        working_dir: str,
        ntasks_per_node: Optional[int] = None,
        gpus_per_node: int = 0
    ) -> List[str]:
        """
        Generate SLURM directives (common for CPU and GPU).
        
        Args:
            num_nodes: Number of nodes
            job_name: Job name
            working_dir: Working directory
            ntasks_per_node: Tasks per node (defaults to cpus_per_node)
            gpus_per_node: GPUs per node (0 for CPU jobs)
        
        Returns:
            List of SLURM directive lines
        """
        # Default ntasks_per_node to full CPU allocation for proper resource accounting
        if ntasks_per_node is None:
            ntasks_per_node = self.config.cpus_per_node
        
        directives = [
            f"#SBATCH --nodes={num_nodes}",
            f"#SBATCH --ntasks-per-node={ntasks_per_node}",
            f"#SBATCH --cpus-per-task=1",
            f"#SBATCH --time={self.config.time_limit}",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH -o {working_dir}/job.out",
            f"#SBATCH -e {working_dir}/job.err",
        ]
        
        if self.config.partition:
            directives.append(f"#SBATCH --partition={self.config.partition}")
        
        if self.config.account:
            directives.append(f"#SBATCH -A {self.config.account}")
        
        if self.config.qos:
            directives.append(f"#SBATCH --qos={self.config.qos}")
        
        if gpus_per_node > 0:
            directives.append(f"#SBATCH --gres=gpu:{gpus_per_node}")
        
        if self.config.exclusive:
            directives.append("#SBATCH --exclusive")
        
        return directives
    
    def generate_module_loads(self) -> str:
        """Generate module load commands."""
        if not self.config.modules:
            return "# No modules specified\n"
        
        lines = ["# Load modules"]
        for module in self.config.modules:
            # Handle both "module load X" and just "X" formats
            if module.startswith("module "):
                lines.append(module)
            else:
                lines.append(f"module load {module}")
        
        return "\n".join(lines) + "\n"
    
    def generate_env_exports(self) -> str:
        """Generate environment variable exports."""
        if not self.config.env_vars:
            return ""
        
        lines = ["# Environment variables"]
        for key, value in self.config.env_vars.items():
            lines.append(f"export {key}={value}")
        
        return "\n".join(lines) + "\n"
    
    def generate_header_comment(
        self,
        num_nodes: int,
        ranks_per_node: int,
        cores_per_rank: int,
        job_type: str = "CPU"
    ) -> str:
        """Generate informative header comment for job script."""
        total_ranks = num_nodes * ranks_per_node
        
        return f'''#!/bin/bash
# =============================================================================
# {job_type} Job Script - Generated by HPC-ScaleTest
# =============================================================================
# Partition: {self.config.partition}
# Job Type: {job_type}
#
# Hardware Configuration (per node):
#   CPU cores: {self.config.cpus_per_node}
#   GPUs: {self.config.gpus_per_node}
#
# MPI Configuration:
#   Total nodes: {num_nodes}
#   MPI ranks per node: {ranks_per_node}
#   Cores per MPI rank: {cores_per_rank}
#   Total MPI ranks: {total_ranks}
#
# Executable: {self.config.executable}
# Input: {self.config.input_file}
# =============================================================================

'''
