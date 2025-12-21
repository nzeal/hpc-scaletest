#!/usr/bin/env python3
"""
MPI Command Generation Module

Generates correct MPI launch commands based on detected topology and
MPI implementation. Supports OpenMPI, Intel MPI, MPICH, and others.

Design Principles:
1. Correct syntax for each MPI implementation
2. Automatic GPU binding via CUDA_VISIBLE_DEVICES
3. No hardcoded assumptions about hardware

This module depends on core/topology.py for hardware detection.

Author: HPC-ScaleTest Contributors
"""

import os
import subprocess
import logging
import re
from dataclasses import dataclass
from typing import Optional, List, Dict, Tuple
from enum import Enum

from core.topology import NodeTopology, MPIMapping, GPUVendor

logger = logging.getLogger(__name__)


class MPIImplementation(Enum):
    """Supported MPI implementations."""
    OPENMPI = "openmpi"
    INTELMPI = "intelmpi"
    MPICH = "mpich"
    MVAPICH = "mvapich"
    CRAY_MPICH = "cray_mpich"
    UNKNOWN = "unknown"


@dataclass
class MPIInfo:
    """Information about the MPI implementation."""
    implementation: MPIImplementation
    version: str
    launcher: str  # mpirun, mpiexec, srun, etc.
    supports_ppr: bool  # --map-by ppr:N:node
    supports_bind_to: bool
    supports_cpus_per_proc: bool
    supports_report_bindings: bool
    
    @classmethod
    def default(cls) -> 'MPIInfo':
        """Return conservative defaults for unknown MPI."""
        return cls(
            implementation=MPIImplementation.UNKNOWN,
            version="",
            launcher="mpirun",
            supports_ppr=False,
            supports_bind_to=False,
            supports_cpus_per_proc=False,
            supports_report_bindings=False,
        )


class MPIDetector:
    """
    Detect MPI implementation and capabilities.
    
    Uses multiple detection methods to identify the MPI implementation
    and determine what command-line options are supported.
    """
    
    def detect(self) -> MPIInfo:
        """
        Detect the current MPI implementation.
        
        Returns:
            MPIInfo with implementation details and capabilities
        """
        # Try mpirun --version
        mpi_info = self._detect_from_mpirun()
        if mpi_info:
            return mpi_info
        
        # Try environment variables
        mpi_info = self._detect_from_env()
        if mpi_info:
            return mpi_info
        
        # Try module system
        mpi_info = self._detect_from_modules()
        if mpi_info:
            return mpi_info
        
        logger.warning("Could not detect MPI implementation, using conservative defaults")
        return MPIInfo.default()
    
    def _detect_from_mpirun(self) -> Optional[MPIInfo]:
        """Detect MPI from mpirun --version output."""
        try:
            result = subprocess.run(
                ['mpirun', '--version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode != 0:
                return None
            
            output = (result.stdout + result.stderr).lower()
            
            # OpenMPI
            if 'open mpi' in output or 'openmpi' in output:
                version_match = re.search(r'(\d+\.\d+\.\d+)', output)
                version = version_match.group(1) if version_match else ""
                return MPIInfo(
                    implementation=MPIImplementation.OPENMPI,
                    version=version,
                    launcher="mpirun",
                    supports_ppr=True,
                    supports_bind_to=True,
                    supports_cpus_per_proc=True,
                    supports_report_bindings=True,
                )
            
            # Intel MPI
            if 'intel' in output:
                version_match = re.search(r'(\d+\.\d+)', output)
                version = version_match.group(1) if version_match else ""
                return MPIInfo(
                    implementation=MPIImplementation.INTELMPI,
                    version=version,
                    launcher="mpirun",
                    supports_ppr=False,  # Intel uses different syntax
                    supports_bind_to=True,
                    supports_cpus_per_proc=False,
                    supports_report_bindings=False,  # Will cause errors
                )
            
            # MPICH
            if 'mpich' in output:
                version_match = re.search(r'(\d+\.\d+)', output)
                version = version_match.group(1) if version_match else ""
                return MPIInfo(
                    implementation=MPIImplementation.MPICH,
                    version=version,
                    launcher="mpiexec",
                    supports_ppr=False,
                    supports_bind_to=False,
                    supports_cpus_per_proc=False,
                    supports_report_bindings=False,
                )
            
            # MVAPICH
            if 'mvapich' in output:
                version_match = re.search(r'(\d+\.\d+)', output)
                version = version_match.group(1) if version_match else ""
                return MPIInfo(
                    implementation=MPIImplementation.MVAPICH,
                    version=version,
                    launcher="mpirun",
                    supports_ppr=False,
                    supports_bind_to=True,
                    supports_cpus_per_proc=False,
                    supports_report_bindings=False,
                )
            
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        return None
    
    def _detect_from_env(self) -> Optional[MPIInfo]:
        """Detect MPI from environment variables."""
        # Check for Intel MPI
        if os.environ.get('I_MPI_ROOT'):
            return MPIInfo(
                implementation=MPIImplementation.INTELMPI,
                version="",
                launcher="mpirun",
                supports_ppr=False,
                supports_bind_to=True,
                supports_cpus_per_proc=False,
                supports_report_bindings=False,
            )
        
        # Check for OpenMPI
        if os.environ.get('OMPI_MCA_PREFIX'):
            return MPIInfo(
                implementation=MPIImplementation.OPENMPI,
                version="",
                launcher="mpirun",
                supports_ppr=True,
                supports_bind_to=True,
                supports_cpus_per_proc=True,
                supports_report_bindings=True,
            )
        
        # Check for Cray MPICH
        if os.environ.get('CRAY_MPICH_DIR'):
            return MPIInfo(
                implementation=MPIImplementation.CRAY_MPICH,
                version="",
                launcher="srun",  # Cray systems typically use srun
                supports_ppr=False,
                supports_bind_to=False,
                supports_cpus_per_proc=False,
                supports_report_bindings=False,
            )
        
        return None
    
    def _detect_from_modules(self) -> Optional[MPIInfo]:
        """Detect MPI from loaded modules."""
        try:
            result = subprocess.run(
                ['bash', '-c', 'module list 2>&1'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            output = result.stdout.lower()
            
            if 'openmpi' in output:
                return MPIInfo(
                    implementation=MPIImplementation.OPENMPI,
                    version="",
                    launcher="mpirun",
                    supports_ppr=True,
                    supports_bind_to=True,
                    supports_cpus_per_proc=True,
                    supports_report_bindings=True,
                )
            
            if 'intelmpi' in output or 'intel-mpi' in output:
                return MPIInfo(
                    implementation=MPIImplementation.INTELMPI,
                    version="",
                    launcher="mpirun",
                    supports_ppr=False,
                    supports_bind_to=True,
                    supports_cpus_per_proc=False,
                    supports_report_bindings=False,
                )
            
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        return None


class MPICommandGenerator:
    """
    Generate MPI launch commands based on topology and MPI implementation.
    
    This class produces the correct mpirun/mpiexec command line for the
    detected MPI implementation and hardware topology.
    
    Usage:
        generator = MPICommandGenerator()
        cmd = generator.generate(topology, mapping, executable="./app")
        # Returns: ["mpirun", "-np", "8", "--map-by", "ppr:4:node", ...]
    """
    
    def __init__(self, mpi_info: Optional[MPIInfo] = None):
        """
        Initialize the command generator.
        
        Args:
            mpi_info: MPI implementation info. If None, auto-detects.
        """
        if mpi_info is None:
            detector = MPIDetector()
            mpi_info = detector.detect()
        self.mpi_info = mpi_info
    
    def generate(
        self,
        topology: NodeTopology,
        mapping: MPIMapping,
        executable: str,
        args: List[str] = None,
        num_nodes: int = 1,
        verbose: bool = False,
        gpu_binding_script: Optional[str] = None,
    ) -> List[str]:
        """
        Generate the MPI launch command.
        
        Args:
            topology: Hardware topology
            mapping: MPI mapping configuration
            executable: Path to executable
            args: Command-line arguments for executable
            num_nodes: Number of nodes
            verbose: Enable verbose/debug output
            gpu_binding_script: Optional GPU binding wrapper script
        
        Returns:
            List of command components (suitable for subprocess)
        """
        args = args or []
        
        # Select generator based on MPI implementation
        if self.mpi_info.implementation == MPIImplementation.OPENMPI:
            return self._generate_openmpi(
                topology, mapping, executable, args, num_nodes, verbose, gpu_binding_script
            )
        elif self.mpi_info.implementation == MPIImplementation.INTELMPI:
            return self._generate_intelmpi(
                topology, mapping, executable, args, num_nodes, verbose, gpu_binding_script
            )
        elif self.mpi_info.implementation == MPIImplementation.CRAY_MPICH:
            return self._generate_cray(
                topology, mapping, executable, args, num_nodes, verbose, gpu_binding_script
            )
        else:
            return self._generate_generic(
                topology, mapping, executable, args, num_nodes, verbose, gpu_binding_script
            )
    
    def _generate_openmpi(
        self,
        topology: NodeTopology,
        mapping: MPIMapping,
        executable: str,
        args: List[str],
        num_nodes: int,
        verbose: bool,
        gpu_binding_script: Optional[str],
    ) -> List[str]:
        """
        Generate OpenMPI command.
        
        OpenMPI syntax with PE (Processing Elements):
            mpirun -np <total_ranks> \
                   --map-by ppr:<ranks_per_node>:node:PE=<cores_per_rank> \
                   [--report-bindings] \
                   [./bind.sh] \
                   <executable> [args]
        
        Example for Leonardo Booster (32 cores, 4 GPUs, 4 nodes):
            mpirun -np 16 --map-by ppr:4:node:PE=8 --report-bindings ./bind.sh <exe>
        """
        cmd = [self.mpi_info.launcher]
        
        # Total number of processes
        total_ranks = mapping.ranks_per_node * num_nodes
        cmd.extend(['-np', str(total_ranks)])
        
        # Process mapping with PE (Processing Elements) for cores per rank
        # Format: ppr:<ranks_per_node>:node:PE=<cores_per_rank>
        cmd.extend(['--map-by', f'ppr:{mapping.ranks_per_node}:node:PE={mapping.cores_per_rank}'])
        
        # Verbose output
        if verbose and self.mpi_info.supports_report_bindings:
            cmd.append('--report-bindings')
        
        # GPU binding wrapper (if provided)
        if gpu_binding_script and topology.gpus > 0:
            cmd.append(gpu_binding_script)
        
        # Executable and arguments
        cmd.append(executable)
        cmd.extend(args)
        
        return cmd
    
    def _generate_intelmpi(
        self,
        topology: NodeTopology,
        mapping: MPIMapping,
        executable: str,
        args: List[str],
        num_nodes: int,
        verbose: bool,
        gpu_binding_script: Optional[str],
    ) -> List[str]:
        """
        Generate Intel MPI command.
        
        Intel MPI syntax:
            mpirun -np <total_ranks> \
                   -ppn <ranks_per_node> \
                   <executable> [args]
        
        Note: Intel MPI does NOT support --report-bindings (will cause failure).
        """
        cmd = [self.mpi_info.launcher]
        
        # Total number of processes
        total_ranks = mapping.ranks_per_node * num_nodes
        cmd.extend(['-np', str(total_ranks)])
        
        # Processes per node
        cmd.extend(['-ppn', str(mapping.ranks_per_node)])
        
        # GPU binding wrapper (if provided)
        if gpu_binding_script and topology.gpus > 0:
            cmd.append(gpu_binding_script)
        
        # Executable and arguments
        cmd.append(executable)
        cmd.extend(args)
        
        return cmd
    
    def _generate_cray(
        self,
        topology: NodeTopology,
        mapping: MPIMapping,
        executable: str,
        args: List[str],
        num_nodes: int,
        verbose: bool,
        gpu_binding_script: Optional[str],
    ) -> List[str]:
        """
        Generate Cray MPICH (srun) command.
        
        Cray systems use srun which gets binding from SLURM:
            srun -n <total_ranks> \
                 --ntasks-per-node=<ranks_per_node> \
                 --cpus-per-task=<cores_per_rank> \
                 <executable> [args]
        """
        cmd = ['srun']
        
        # Total number of tasks
        total_ranks = mapping.ranks_per_node * num_nodes
        cmd.extend(['-n', str(total_ranks)])
        
        # Tasks per node
        cmd.extend(['--ntasks-per-node', str(mapping.ranks_per_node)])
        
        # CPUs per task
        if mapping.cores_per_rank > 1:
            cmd.extend(['--cpus-per-task', str(mapping.cores_per_rank)])
        
        # GPU binding (srun handles this via --gpus-per-task if needed)
        if topology.gpus > 0 and mapping.gpus_per_rank > 0:
            cmd.extend(['--gpus-per-task', str(mapping.gpus_per_rank)])
        
        # Executable and arguments
        cmd.append(executable)
        cmd.extend(args)
        
        return cmd
    
    def _generate_generic(
        self,
        topology: NodeTopology,
        mapping: MPIMapping,
        executable: str,
        args: List[str],
        num_nodes: int,
        verbose: bool,
        gpu_binding_script: Optional[str],
    ) -> List[str]:
        """
        Generate generic MPI command (minimal options).
        
        Uses only widely-supported options:
            mpirun -np <total_ranks> <executable> [args]
        """
        cmd = [self.mpi_info.launcher]
        
        # Total number of processes
        total_ranks = mapping.ranks_per_node * num_nodes
        cmd.extend(['-np', str(total_ranks)])
        
        # GPU binding wrapper (if provided)
        if gpu_binding_script and topology.gpus > 0:
            cmd.append(gpu_binding_script)
        
        # Executable and arguments
        cmd.append(executable)
        cmd.extend(args)
        
        return cmd
    
    def generate_command_string(
        self,
        topology: NodeTopology,
        mapping: MPIMapping,
        executable: str,
        args: List[str] = None,
        num_nodes: int = 1,
        verbose: bool = False,
        gpu_binding_script: Optional[str] = None,
    ) -> str:
        """
        Generate command as a shell string.
        
        Convenience method that returns a single string suitable for
        inclusion in shell scripts.
        """
        cmd = self.generate(
            topology, mapping, executable, args, num_nodes, verbose, gpu_binding_script
        )
        
        # Quote arguments that need it
        def quote_if_needed(s: str) -> str:
            if ' ' in s or '"' in s or "'" in s or '$' in s:
                return f'"{s}"'
            return s
        
        return ' '.join(quote_if_needed(c) for c in cmd)


def generate_gpu_binding_script(
    topology: NodeTopology,
    output_dir: str = ".",
    script_name: str = "bind_gpu.sh"
) -> str:
    """
    Generate a GPU binding wrapper script.
    
    This script uses OMPI_COMM_WORLD_LOCAL_RANK (OpenMPI) or
    PMI_LOCAL_RANK (MPICH/Intel) to set CUDA_VISIBLE_DEVICES
    such that each MPI rank gets a unique GPU.
    
    Args:
        topology: Hardware topology (for GPU vendor)
        output_dir: Directory to write script
        script_name: Name of the script file
    
    Returns:
        Path to the generated script
    """
    # Determine the environment variable for GPU visibility
    if topology.gpu_vendor == GPUVendor.NVIDIA:
        gpu_env_var = "CUDA_VISIBLE_DEVICES"
    elif topology.gpu_vendor == GPUVendor.AMD:
        gpu_env_var = "ROCR_VISIBLE_DEVICES"
    elif topology.gpu_vendor == GPUVendor.INTEL:
        gpu_env_var = "ZE_AFFINITY_MASK"
    else:
        gpu_env_var = "CUDA_VISIBLE_DEVICES"  # Default to NVIDIA
    
    script_content = f'''#!/bin/bash
# GPU Binding Script
# Generated by HPC-ScaleTest
#
# This script binds each MPI rank to a unique GPU using the local rank.
# It supports OpenMPI (OMPI_COMM_WORLD_LOCAL_RANK) and 
# MPICH/Intel MPI (PMI_LOCAL_RANK, MPI_LOCALRANKID).
#
# Usage: mpirun ... ./bind_gpu.sh <executable> [args]

# Determine local rank from available environment variables
if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
    LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ -n "$PMI_LOCAL_RANK" ]; then
    LOCAL_RANK=$PMI_LOCAL_RANK
elif [ -n "$MPI_LOCALRANKID" ]; then
    LOCAL_RANK=$MPI_LOCALRANKID
elif [ -n "$SLURM_LOCALID" ]; then
    LOCAL_RANK=$SLURM_LOCALID
else
    # Fallback: cannot determine local rank
    echo "WARNING: Cannot determine local MPI rank for GPU binding" >&2
    LOCAL_RANK=0
fi

# Set GPU visibility to only the GPU for this rank
export {gpu_env_var}=$LOCAL_RANK

# Debug output (uncomment if needed)
# echo "Rank ${{OMPI_COMM_WORLD_RANK:-$PMI_RANK}}: LOCAL_RANK=$LOCAL_RANK, {gpu_env_var}=${gpu_env_var}"

# Execute the actual command
exec "$@"
'''
    
    script_path = os.path.join(output_dir, script_name)
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    # Make executable
    os.chmod(script_path, 0o755)
    
    logger.info(f"Generated GPU binding script: {script_path}")
    return script_path


def get_mpi_command(
    topology: NodeTopology,
    mapping: MPIMapping,
    executable: str,
    args: List[str] = None,
    num_nodes: int = 1,
    verbose: bool = False,
    include_gpu_binding: bool = True,
    output_dir: str = ".",
) -> Tuple[List[str], Optional[str]]:
    """
    Convenience function to generate MPI command with optional GPU binding.
    
    Args:
        topology: Hardware topology
        mapping: MPI mapping configuration
        executable: Path to executable
        args: Command-line arguments
        num_nodes: Number of nodes
        verbose: Enable verbose output
        include_gpu_binding: Generate GPU binding script if GPUs present
        output_dir: Directory for generated scripts
    
    Returns:
        Tuple of (command_list, gpu_binding_script_path)
    """
    gpu_binding_script = None
    
    # Generate GPU binding script if needed
    if include_gpu_binding and topology.gpus > 0:
        gpu_binding_script = generate_gpu_binding_script(topology, output_dir)
    
    # Generate MPI command
    generator = MPICommandGenerator()
    cmd = generator.generate(
        topology=topology,
        mapping=mapping,
        executable=executable,
        args=args,
        num_nodes=num_nodes,
        verbose=verbose,
        gpu_binding_script=gpu_binding_script,
    )
    
    return cmd, gpu_binding_script
