"""
MPI Run Launcher Backend

This launcher uses the centralized topology detection (core/topology.py) and
MPI command generation (core/mpi_command.py) modules to produce correct
mpirun commands without hardcoded assumptions.

Design Principles:
1. All topology detection delegated to core/topology.py
2. All MPI command generation delegated to core/mpi_command.py
3. No hardcoded values in this module
4. Graceful fallback if new modules not available

Author: HPC-ScaleTest Contributors
"""

import logging
import os
from typing import List, Optional
from pathlib import Path

from core.abstracts import LauncherInterface
from core.config import JobConfig, ResourceConfig

logger = logging.getLogger(__name__)

# Import centralized topology and MPI command modules
try:
    from core.topology import (
        TopologyDetector, NodeTopology, MPIMapping,
        get_topology_detector, GPUVendor
    )
    from core.mpi_command import (
        MPICommandGenerator, MPIDetector, MPIInfo,
        generate_gpu_binding_script, MPIImplementation
    )
    HAS_TOPOLOGY_MODULE = True
except ImportError:
    HAS_TOPOLOGY_MODULE = False
    logger.warning("Centralized topology module not available")

# Legacy imports for fallback
try:
    from utils.mpi_detector import MPIDetector as LegacyMPIDetector
    HAS_LEGACY_MPI_DETECTOR = True
except ImportError:
    HAS_LEGACY_MPI_DETECTOR = False

try:
    from utils.gpu_bind_generator import generate_gpu_bind_script
    HAS_LEGACY_GPU_BIND = True
except ImportError:
    HAS_LEGACY_GPU_BIND = False


class MpiRunLauncher(LauncherInterface):
    """
    MPI Launcher with automatic topology detection and correct command generation.
    
    This launcher automatically:
    1. Detects hardware topology (CPUs, GPUs per node)
    2. Computes optimal MPI mapping (ranks, cores per rank)
    3. Generates correct mpirun command with proper syntax
    4. Handles GPU binding via CUDA_VISIBLE_DEVICES
    
    The launcher supports:
    - OpenMPI: --map-by ppr:N:node --bind-to core --cpus-per-proc M
    - Intel MPI: -ppn N
    - MPICH: -np N
    - Generic fallback for unknown implementations
    
    Usage:
        launcher = MpiRunLauncher(options={'verbose': True})
        cmd = launcher.generate_launch_command(job_config, executable, resource_config)
    """
    
    def __init__(self, options: dict = None):
        """
        Initialize the MPI launcher.
        
        Args:
            options: Optional configuration dict with keys:
                - launcher: Override launcher name (default: auto-detect)
                - verbose: Enable verbose output (default: False)
                - gpu_binding: Enable GPU binding script (default: True for GPU jobs)
        """
        super().__init__(options)
        self._topology_detector = None
        self._mpi_detector = None
        self._mpi_info = None
    
    @property
    def topology_detector(self) -> Optional['TopologyDetector']:
        """Get or create topology detector."""
        if self._topology_detector is None and HAS_TOPOLOGY_MODULE:
            self._topology_detector = get_topology_detector()
        return self._topology_detector
    
    @property
    def mpi_info(self) -> Optional['MPIInfo']:
        """Get or detect MPI information."""
        if self._mpi_info is None and HAS_TOPOLOGY_MODULE:
            self._mpi_detector = MPIDetector()
            self._mpi_info = self._mpi_detector.detect()
        return self._mpi_info
    
    def supports_gpu_binding(self) -> bool:
        """Check if launcher supports GPU binding."""
        return True
    
    def generate_launch_command(
        self,
        job_config: JobConfig,
        executable: List[str],
        resource_config: ResourceConfig
    ) -> List[str]:
        """
        Generate mpirun launch command with automatic topology detection.
        
        This method:
        1. Detects hardware topology for the partition
        2. Computes optimal MPI mapping
        3. Generates correct mpirun command for detected MPI implementation
        4. Includes GPU binding script if GPUs are present
        
        Args:
            job_config: Job configuration
            executable: Executable and arguments as list
            resource_config: Resource configuration
        
        Returns:
            List of command components (launcher + arguments + executable)
        """
        if HAS_TOPOLOGY_MODULE:
            return self._generate_with_topology(job_config, executable, resource_config)
        else:
            return self._generate_legacy(job_config, executable, resource_config)
    
    def _generate_with_topology(
        self,
        job_config: JobConfig,
        executable: List[str],
        resource_config: ResourceConfig
    ) -> List[str]:
        """Generate command using centralized topology detection."""
        
        # Step 1: Detect or use provided topology
        partition = getattr(resource_config, 'partition', None)
        
        try:
            topology = self.topology_detector.detect(partition)
        except Exception as e:
            logger.warning(f"Topology detection failed: {e}, using resource_config values")
            # Create topology from resource_config
            topology = NodeTopology(
                cpu_cores=resource_config.procs_per_node,
                gpus=getattr(resource_config, 'gpus_per_node', 0),
                gpu_vendor=GPUVendor.NVIDIA if getattr(resource_config, 'gpus_per_node', 0) > 0 else GPUVendor.NONE,
            )
        
        # Step 2: Determine MPI mapping
        # Check for user overrides in resource_config
        user_ranks = getattr(resource_config, 'actual_mpi_tasks', None)
        user_cores = getattr(resource_config, 'cores_per_task', None)
        
        try:
            mapping = self.topology_detector.compute_mpi_mapping(
                topology=topology,
                num_nodes=job_config.num_nodes,
                user_ranks_per_node=user_ranks,
                user_cores_per_rank=user_cores,
            )
        except ValueError as e:
            logger.error(f"Invalid MPI mapping: {e}")
            raise
        
        # Step 3: Generate MPI command
        mpi_info = self.mpi_info
        
        # Override launcher if specified
        launcher = self.options.get('launcher')
        if launcher:
            mpi_info = MPIInfo(
                implementation=mpi_info.implementation if mpi_info else MPIImplementation.UNKNOWN,
                version=mpi_info.version if mpi_info else "",
                launcher=launcher,
                supports_ppr=mpi_info.supports_ppr if mpi_info else False,
                supports_bind_to=mpi_info.supports_bind_to if mpi_info else False,
                supports_cpus_per_proc=mpi_info.supports_cpus_per_proc if mpi_info else False,
                supports_report_bindings=mpi_info.supports_report_bindings if mpi_info else False,
            )
        
        generator = MPICommandGenerator(mpi_info)
        
        # Determine if we need GPU binding
        gpu_binding_script = None
        enable_gpu_binding = self.options.get('gpu_binding', True)
        
        if enable_gpu_binding and topology.gpus > 0:
            # Generate GPU binding script inline
            gpu_binding_script = './bind.sh'
        
        verbose = self.options.get('verbose', False)
        
        # Generate command
        exe_path = executable[0] if executable else ''
        exe_args = executable[1:] if len(executable) > 1 else []
        
        cmd = generator.generate(
            topology=topology,
            mapping=mapping,
            executable=exe_path,
            args=exe_args,
            num_nodes=job_config.num_nodes,
            verbose=verbose,
            gpu_binding_script=gpu_binding_script,
        )
        
        # Log the command
        logger.info(f"Generated MPI command for {mpi_info.implementation.value if mpi_info else 'unknown'} MPI:")
        logger.info(f"  Topology: {topology.cpu_cores} CPUs, {topology.gpus} GPUs per node")
        logger.info(f"  Mapping: {mapping.ranks_per_node} ranks/node × {mapping.cores_per_rank} cores/rank")
        logger.info(f"  Command: {' '.join(cmd)}")
        
        return cmd
    
    def _generate_legacy(
        self,
        job_config: JobConfig,
        executable: List[str],
        resource_config: ResourceConfig
    ) -> List[str]:
        """
        Legacy command generation (fallback if topology module unavailable).
        
        This provides backward compatibility but may not generate optimal commands.
        """
        logger.warning("Using legacy MPI command generation (topology module not available)")
        
        launcher = self.options.get('launcher', 'mpirun')
        cmd = [launcher]
        
        # Determine task count
        if hasattr(resource_config, 'actual_mpi_tasks') and resource_config.actual_mpi_tasks:
            tasks_per_node = resource_config.actual_mpi_tasks
            total_tasks = tasks_per_node * job_config.num_nodes
            cores_per_task = getattr(resource_config, 'cores_per_task', 1)
            is_gpu_job = True
        else:
            total_tasks = job_config.num_nodes * resource_config.procs_per_node
            tasks_per_node = resource_config.procs_per_node
            cores_per_task = 1
            is_gpu_job = getattr(resource_config, 'gpus_per_node', 0) > 0
        
        # Number of processes
        cmd.extend(['-np', str(total_tasks)])
        
        # Process mapping for GPU jobs with PE (Processing Elements)
        if is_gpu_job:
            # Use ppr:X:node:PE=Y syntax for process mapping with cores per rank
            cmd.extend(['--map-by', f'ppr:{tasks_per_node}:node:PE={cores_per_task}'])
        
        # Add executable
        cmd.extend(executable)
        
        return cmd
    
    def get_gpu_binding_script_content(self, gpu_vendor: str = 'nvidia') -> str:
        """
        Get the content of the GPU binding script.
        
        This can be used by job generators that need to embed the script
        directly in the job script.
        
        Args:
            gpu_vendor: GPU vendor ('nvidia', 'amd', 'intel')
        
        Returns:
            Shell script content as string
        """
        env_vars = {
            'nvidia': 'CUDA_VISIBLE_DEVICES',
            'amd': 'ROCR_VISIBLE_DEVICES',
            'intel': 'ZE_AFFINITY_MASK',
        }
        gpu_env = env_vars.get(gpu_vendor.lower(), 'CUDA_VISIBLE_DEVICES')
        
        return f'''#!/bin/bash
# GPU Binding Script - Generated by HPC-ScaleTest
# Binds each MPI rank to a unique GPU using local rank

if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
    LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ -n "$PMI_LOCAL_RANK" ]; then
    LOCAL_RANK=$PMI_LOCAL_RANK
elif [ -n "$MPI_LOCALRANKID" ]; then
    LOCAL_RANK=$MPI_LOCALRANKID
elif [ -n "$SLURM_LOCALID" ]; then
    LOCAL_RANK=$SLURM_LOCALID
else
    LOCAL_RANK=0
fi

export {gpu_env}=$LOCAL_RANK
exec "$@"
'''


# Backward compatibility alias
MpirunLauncher = MpiRunLauncher
