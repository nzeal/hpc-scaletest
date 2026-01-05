"""
MPI Run Launcher Backend

This launcher uses the unified execution module (core/unified_execution.py)
to generate correct mpirun commands with automatic topology detection.

Key Features:
=============
- Automatic topology detection (CPUs, GPUs per node)
- Correct mpirun syntax: -np N --map-by ppr:X:node:PE=Y
- GPU binding via CUDA_VISIBLE_DEVICES (NVIDIA) or ROCR_VISIBLE_DEVICES (AMD)
- Support for OpenMPI, Intel MPI, MPICH, and other implementations

Design Principles:
==================
- NO HARDCODED VALUES - All topology detected at runtime
- SYSTEM AGNOSTIC - Works on any HPC system
- SINGLE SOURCE OF TRUTH - Uses core/unified_execution.py

Author: HPC-ScaleTest Contributors
"""

import logging
import os
from typing import List, Optional

from core.abstracts import LauncherInterface
from core.config import JobConfig, ResourceConfig

logger = logging.getLogger(__name__)

# Import unified execution module
try:
    from core.unified_execution import (
        UnifiedExecutor, TopologyDetector, HardwareTopology, GPUVendor
    )
    HAS_UNIFIED_EXECUTION = True
except ImportError:
    HAS_UNIFIED_EXECUTION = False
    logger.warning("Unified execution module not available")


class MpiRunLauncher(LauncherInterface):
    """
    MPI Launcher with automatic topology detection and correct command generation.
    
    This launcher automatically:
    1. Detects hardware topology (CPUs, GPUs per node)
    2. Computes optimal MPI mapping (ranks, cores per rank)
    3. Generates correct mpirun command
    4. Handles GPU binding via bind.sh
    
    MPI Command Format (OpenMPI):
    =============================
    mpirun -np <total_ranks> --map-by ppr:<ranks_per_node>:node:PE=<cores_per_rank> [--report-bindings] [./bind.sh] <executable> [args]
    
    Example (Leonardo Booster, 4 nodes):
    - 32 CPU cores per node, 4 GPUs per node
    - ranks_per_node = 4 (1 per GPU)
    - cores_per_rank = 8 (32/4)
    - Command: mpirun -np 16 --map-by ppr:4:node:PE=8 --report-bindings ./bind.sh $BINARY/iPIC3D os-stdin
    """
    
    def __init__(self, options: dict = None):
        """
        Initialize the MPI launcher.
        
        Args:
            options: Optional configuration dict with keys:
                - launcher: Override launcher name (default: mpirun)
                - verbose: Enable verbose output (default: False)
                - gpu_binding: Enable GPU binding script (default: True for GPU jobs)
                - report_bindings: Include --report-bindings (default: True)
        """
        super().__init__(options)
        self._executor: Optional[UnifiedExecutor] = None
    
    def _get_executor(self, partition: str = "") -> UnifiedExecutor:
        """Get or create unified executor."""
        if not HAS_UNIFIED_EXECUTION:
            raise RuntimeError("Unified execution module not available")
        
        if self._executor is None or (partition and self._executor.partition != partition):
            self._executor = UnifiedExecutor(partition=partition)
        return self._executor
    
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
        2. Computes MPI mapping (ranks per node, cores per rank)
        3. Generates correct mpirun command
        
        Args:
            job_config: Job configuration
            executable: Executable and arguments as list
            resource_config: Resource configuration
        
        Returns:
            List of command components (launcher + arguments + executable)
        
        Example output for Leonardo (32 cores, 4 GPUs), 4 nodes:
            ['mpirun', '-np', '16', '--map-by', 'ppr:4:node:PE=8', 
             '--report-bindings', './bind.sh', '$BINARY/iPIC3D', 'os-stdin']
        """
        if HAS_UNIFIED_EXECUTION:
            return self._generate_with_unified_execution(
                job_config, executable, resource_config
            )
        else:
            return self._generate_legacy(job_config, executable, resource_config)
    
    def _generate_with_unified_execution(
        self,
        job_config: JobConfig,
        executable: List[str],
        resource_config: ResourceConfig
    ) -> List[str]:
        """Generate command using unified execution module."""
        
        # Get partition
        partition = getattr(resource_config, 'partition', None) or ''
        
        # Get executor and detect topology
        executor = self._get_executor(partition)
        topology = executor.detect_topology()
        
        # Check for user overrides
        user_gpus = getattr(resource_config, 'gpus_per_node', 0)
        if user_gpus > topology.gpus_per_node:
            # User configured more GPUs - rebuild topology
            topology = HardwareTopology(
                cpu_cores_per_node=topology.cpu_cores_per_node,
                gpus_per_node=user_gpus,
                gpu_vendor=topology.gpu_vendor if topology.gpus_per_node > 0 else GPUVendor.NVIDIA,
                partition=partition,
                detection_method=topology.detection_method + " + user_override"
            )
            executor._topology = topology
        
        # Calculate MPI layout
        total_ranks = job_config.num_nodes * topology.ranks_per_node
        
        # Build command
        launcher = self.options.get('launcher', 'mpirun')
        report_bindings = self.options.get('report_bindings', True)
        gpu_binding = self.options.get('gpu_binding', True)
        
        cmd = [
            launcher,
            '-np', str(total_ranks),
            '--map-by', f'ppr:{topology.ranks_per_node}:node:PE={topology.cores_per_rank}',
            '--report-bindings'
        ]
        
        # GPU binding wrapper - use bind.sh for GPU jobs
        if gpu_binding and topology.gpus_per_node > 0:
            cmd.append('./bind.sh')
        
        # Add executable and arguments
        cmd.extend(executable)
        
        logger.info(f"Generated MPI command:")
        logger.info(f"  Topology: {topology.cpu_cores_per_node} CPUs, {topology.gpus_per_node} GPUs per node")
        logger.info(f"  Mapping: {topology.ranks_per_node} ranks/node × {topology.cores_per_rank} cores/rank")
        logger.info(f"  Total ranks: {total_ranks}")
        logger.info(f"  Command: {' '.join(cmd)}")
        
        return cmd
    
    def _generate_legacy(
        self,
        job_config: JobConfig,
        executable: List[str],
        resource_config: ResourceConfig
    ) -> List[str]:
        """
        Legacy command generation (fallback if unified execution unavailable).
        
        This uses resource_config values directly without automatic detection.
        """
        logger.warning("Using legacy MPI command generation")
        
        launcher = self.options.get('launcher', 'mpirun')
        cmd = [launcher]
        
        # Determine MPI layout from resource config
        gpus_per_node = getattr(resource_config, 'gpus_per_node', 0)
        procs_per_node = resource_config.procs_per_node
        
        if gpus_per_node > 0:
            # GPU job: 1 rank per GPU
            ranks_per_node = gpus_per_node
            cores_per_rank = procs_per_node // gpus_per_node if gpus_per_node > 0 else 1
            total_ranks = job_config.num_nodes * ranks_per_node
            
            cmd.extend(['-np', str(total_ranks)])
            cmd.extend(['--map-by', f'ppr:{ranks_per_node}:node:PE={cores_per_rank}'])
            
            if self.options.get('report_bindings', True):
                cmd.append('--report-bindings')
            
            if self.options.get('gpu_binding', True):
                cmd.append('./bind.sh')
        else:
            # CPU job: all processes
            total_procs = job_config.num_nodes * procs_per_node
            cmd.extend(['-np', str(total_procs)])
        
        # Add executable
        cmd.extend(executable)
        
        return cmd
    
    def get_gpu_binding_script_content(self, gpu_vendor: str = 'nvidia') -> str:
        """
        Get the content of the GPU binding script.
        
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
elif [ -n "$MPI_LOCALRANKID" ]; then
    LOCAL_RANK=$MPI_LOCALRANKID
elif [ -n "$PMI_LOCAL_RANK" ]; then
    LOCAL_RANK=$PMI_LOCAL_RANK
elif [ -n "$SLURM_LOCALID" ]; then
    LOCAL_RANK=$SLURM_LOCALID
elif [ -n "$MV2_COMM_WORLD_LOCAL_RANK" ]; then
    LOCAL_RANK=$MV2_COMM_WORLD_LOCAL_RANK
else
    LOCAL_RANK=0
fi

export {gpu_env}=$LOCAL_RANK
exec "$@"
'''


# Backward compatibility alias
MpirunLauncher = MpiRunLauncher
