#!/usr/bin/env python3
"""
Comprehensive Hardware Detection Module

Automatically detects system hardware and computes optimal MPI task mapping
for GPU-accelerated HPC applications.

Features:
- Auto-detects GPUs (NVIDIA, AMD, Intel)
- Auto-detects CPU cores per node
- Computes optimal --map-by ppr:X:node --bind-to core --cpus-per-proc Y mapping
- System-agnostic (works on Leonardo, LUMI, and other systems)
"""

import subprocess
import logging
import os
from typing import Optional, Tuple, Dict
from dataclasses import dataclass

# Import GPU detector
try:
    from utils.gpu_detection import GPUDetector, GPUInfo
    HAS_GPU_DETECTOR = True
except ImportError:
    HAS_GPU_DETECTOR = False
    GPUInfo = None

logger = logging.getLogger(__name__)


@dataclass
class HardwareConfig:
    """Complete hardware configuration for a compute node."""
    # CPU information
    cores_per_node: int
    threads_per_core: int
    physical_cores_per_node: int
    
    # GPU information
    gpus_per_node: int
    gpu_vendor: Optional[str] = None
    gpu_model: Optional[str] = None
    
    # Computed MPI mapping
    mpi_tasks_per_node: int = 0  # X in --map-by ppr:X:node
    cores_per_task: int = 0      # Y for --cpus-per-proc Y
    
    def get_ppr_mapping(self) -> str:
        """
        Generate the --map-by ppr mapping string.
        
        Returns:
            String like "ppr:4:node" for use with mpirun --map-by
        """
        if self.mpi_tasks_per_node > 0:
            return f"ppr:{self.mpi_tasks_per_node}:node"
        return ""
    
    def get_mpirun_args(self) -> list:
        """
        Generate complete mpirun mapping arguments.
        
        Returns:
            List of arguments for mpirun command
        """
        args = []
        if self.mpi_tasks_per_node > 0:
            args.extend(["--map-by", f"ppr:{self.mpi_tasks_per_node}:node"])
            args.extend(["--bind-to", "core"])
            if self.cores_per_task > 1:
                args.extend(["--cpus-per-proc", str(self.cores_per_task)])
        return args
    
    def __str__(self) -> str:
        """Human-readable hardware summary."""
        lines = [
            f"Hardware Configuration:",
            f"  CPU: {self.cores_per_node} cores/node ({self.physical_cores_per_node} physical)",
            f"  GPU: {self.gpus_per_node} {self.gpu_vendor or 'N/A'} GPUs/node"
        ]
        
        if self.gpus_per_node > 0:
            lines.append(f"  Model: {self.gpu_model}")
            lines.append(f"\nOptimal MPI Mapping:")
            lines.append(f"  Tasks per node: {self.mpi_tasks_per_node} (1 per GPU)")
            lines.append(f"  Cores per task: {self.cores_per_task}")
            mpi_args = self.get_mpirun_args()
            lines.append(f"  Command: mpirun {' '.join(mpi_args)}")
        
        return "\n".join(lines)


class HardwareDetector:
    """
    Automatically detect system hardware and compute optimal configuration.
    
    This class combines CPU and GPU detection to automatically determine
    the correct MPI task mapping for GPU-accelerated applications.
    """
    
    def __init__(self):
        self.gpu_detector = GPUDetector() if HAS_GPU_DETECTOR else None
        self.hardware_config: Optional[HardwareConfig] = None
    
    def detect_hardware(self, partition_name: Optional[str] = None) -> HardwareConfig:
        """
        Automatically detect all hardware and compute optimal configuration.
        
        Args:
            partition_name: Optional SLURM partition name for partition-specific detection
            
        Returns:
            HardwareConfig with complete hardware info and computed MPI mapping
        """
        # Detect CPUs
        cores_per_node, threads_per_core, physical_cores = self._detect_cpu_cores()
        
        logger.info(f"✓ Detected CPU hardware:")
        logger.info(f"  Cores per node: {cores_per_node}")
        logger.info(f"  Threads per core: {threads_per_core}")
        logger.info(f"  Physical cores: {physical_cores}")
        
        # Detect GPUs
        gpu_info = None
        gpus_per_node = 0
        gpu_vendor = None
        gpu_model = None
        
        if self.gpu_detector:
            gpu_info = self.gpu_detector.detect(partition_name=partition_name)
            if gpu_info:
                gpus_per_node = gpu_info.count_per_node
                gpu_vendor = gpu_info.vendor
                gpu_model = gpu_info.model
        
        # Compute optimal MPI mapping
        if gpus_per_node > 0:
            # GPU system: 1 MPI task per GPU
            mpi_tasks_per_node = gpus_per_node
            cores_per_task = cores_per_node // gpus_per_node
            
            logger.info(f"\n✓ Computed GPU-aware MPI mapping:")
            logger.info(f"  MPI tasks per node: {mpi_tasks_per_node} (1 per GPU)")
            logger.info(f"  Cores per task: {cores_per_task}")
            logger.info(f"  Mapping: --map-by ppr:{mpi_tasks_per_node}:node --bind-to core --cpus-per-proc {cores_per_task}")
        else:
            # CPU-only system: use all cores
            mpi_tasks_per_node = cores_per_node
            cores_per_task = 1
            
            logger.info(f"\n✓ CPU-only system (no GPUs detected)")
            logger.info(f"  MPI tasks per node: {mpi_tasks_per_node}")
        
        # Create hardware config
        config = HardwareConfig(
            cores_per_node=cores_per_node,
            threads_per_core=threads_per_core,
            physical_cores_per_node=physical_cores,
            gpus_per_node=gpus_per_node,
            gpu_vendor=gpu_vendor,
            gpu_model=gpu_model,
            mpi_tasks_per_node=mpi_tasks_per_node,
            cores_per_task=cores_per_task
        )
        
        self.hardware_config = config
        return config
    
    def _detect_cpu_cores(self) -> Tuple[int, int, int]:
        """
        Detect CPU cores per node using multiple methods.
        
        Returns:
            Tuple of (total_cores, threads_per_core, physical_cores)
        """
        # Method 1: Try lscpu (most reliable on Linux)
        cores, threads, physical = self._detect_with_lscpu()
        if cores > 0:
            logger.debug(f"lscpu: {cores} cores, {threads} threads/core")
            return cores, threads, physical
        
        # Method 2: Try /proc/cpuinfo
        cores, threads, physical = self._detect_from_cpuinfo()
        if cores > 0:
            logger.debug(f"cpuinfo: {cores} cores")
            return cores, threads, physical
        
        # Method 3: Try SLURM environment variables
        cores, threads, physical = self._detect_from_slurm_env()
        if cores > 0:
            logger.debug(f"SLURM env: {cores} cores")
            return cores, threads, physical
        
        # Method 4: Fallback to os.cpu_count()
        cores = os.cpu_count() or 1
        logger.warning(f"Using os.cpu_count() fallback: {cores} cores")
        return cores, 1, cores
    
    def _detect_with_lscpu(self) -> Tuple[int, int, int]:
        """
        Detect CPU info using lscpu command.
        
        Returns:
            (total_cores, threads_per_core, physical_cores)
        """
        try:
            result = subprocess.run(
                ['lscpu'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode != 0:
                return 0, 0, 0
            
            cores_per_socket = 1
            sockets = 1
            threads_per_core = 1
            
            for line in result.stdout.split('\n'):
                line = line.strip()
                
                if line.startswith('CPU(s):'):
                    # Total logical CPUs
                    total_cpus = int(line.split(':')[1].strip())
                elif line.startswith('Thread(s) per core:'):
                    threads_per_core = int(line.split(':')[1].strip())
                elif line.startswith('Core(s) per socket:'):
                    cores_per_socket = int(line.split(':')[1].strip())
                elif line.startswith('Socket(s):'):
                    sockets = int(line.split(':')[1].strip())
            
            # Calculate values
            physical_cores = cores_per_socket * sockets
            total_cores = physical_cores * threads_per_core
            
            return total_cores, threads_per_core, physical_cores
            
        except FileNotFoundError:
            logger.debug("lscpu command not found")
            return 0, 0, 0
        except Exception as e:
            logger.debug(f"lscpu detection failed: {e}")
            return 0, 0, 0
    
    def _detect_from_cpuinfo(self) -> Tuple[int, int, int]:
        """
        Detect CPU info from /proc/cpuinfo.
        
        Returns:
            (total_cores, threads_per_core, physical_cores)
        """
        try:
            with open('/proc/cpuinfo', 'r') as f:
                content = f.read()
            
            # Count processor entries
            processor_count = content.count('processor\t:')
            
            # Try to detect siblings (hyperthreading)
            siblings = 1
            cores = 1
            for line in content.split('\n'):
                if line.startswith('siblings'):
                    siblings = int(line.split(':')[1].strip())
                elif line.startswith('cpu cores'):
                    cores = int(line.split(':')[1].strip())
            
            threads_per_core = siblings // cores if cores > 0 else 1
            physical_cores = processor_count // threads_per_core
            
            return processor_count, threads_per_core, physical_cores
            
        except Exception as e:
            logger.debug(f"cpuinfo detection failed: {e}")
            return 0, 0, 0
    
    def _detect_from_slurm_env(self) -> Tuple[int, int, int]:
        """
        Detect CPU info from SLURM environment variables.
        
        Returns:
            (total_cores, threads_per_core, physical_cores)
        """
        try:
            # SLURM sets these when job is running
            if 'SLURM_CPUS_ON_NODE' in os.environ:
                cores = int(os.environ['SLURM_CPUS_ON_NODE'])
                return cores, 1, cores  # Assume no hyperthreading info
            
            if 'SLURM_JOB_CPUS_PER_NODE' in os.environ:
                # Format can be like "32" or "32(x2)" for 2 nodes
                cpus_str = os.environ['SLURM_JOB_CPUS_PER_NODE']
                cores = int(cpus_str.split('(')[0])
                return cores, 1, cores
            
        except Exception as e:
            logger.debug(f"SLURM env detection failed: {e}")
        
        return 0, 0, 0
    
    def configure_resource_config(self, resource_config):
        """
        Automatically configure a ResourceConfig object with detected hardware.
        
        IMPORTANT: This method PRESERVES user-configured values.
        If gpus_per_node was already set by the user, it will NOT be overwritten.
        
        Args:
            resource_config: ResourceConfig instance to configure
        """
        if not self.hardware_config:
            self.detect_hardware()
        
        config = self.hardware_config
        
        # CRITICAL FIX: Preserve user's gpus_per_node if already set
        # Only use auto-detected value if user didn't specify one
        user_gpus = getattr(resource_config, 'gpus_per_node', 0) or 0
        if user_gpus > 0:
            # User already set gpus_per_node - preserve it
            logger.info(f"  Preserving user-configured gpus_per_node: {user_gpus}")
            effective_gpus = user_gpus
        else:
            # Use auto-detected value
            resource_config.gpus_per_node = config.gpus_per_node
            effective_gpus = config.gpus_per_node
        
        # Configure GPU tasks if GPUs are present
        if effective_gpus > 0:
            # Calculate MPI tasks based on effective GPU count
            # If user preserved their gpus_per_node, compute tasks from that
            if user_gpus > 0 and config.gpus_per_node == 0:
                # User specified GPUs but auto-detect found none
                # Use user's cpus_per_node if available, otherwise use auto-detected or default
                user_cpus = getattr(resource_config, 'cpus_per_node', 0) or 0
                if user_cpus > 0:
                    cores_per_node = user_cpus
                elif config.cores_per_node > 0:
                    cores_per_node = config.cores_per_node
                else:
                    cores_per_node = 32  # Default for modern GPU nodes
                
                mpi_tasks = effective_gpus
                cores_per_task = cores_per_node // effective_gpus
                logger.info(f"  Computing MPI layout from user config:")
                logger.info(f"    Using {cores_per_node} cores / {effective_gpus} GPUs = {cores_per_task} cores/task")
            else:
                # Use auto-detected values
                mpi_tasks = config.mpi_tasks_per_node
                cores_per_task = config.cores_per_task
            
            resource_config.actual_mpi_tasks = mpi_tasks
            resource_config.cores_per_task = cores_per_task
            
            logger.info(f"✓ Auto-configured ResourceConfig for GPU execution:")
            logger.info(f"  gpus_per_node: {effective_gpus}")
            logger.info(f"  actual_mpi_tasks: {resource_config.actual_mpi_tasks}")
            logger.info(f"  cores_per_task: {resource_config.cores_per_task}")
            logger.info(f"  → MPI mapping: --map-by ppr:{mpi_tasks}:node --bind-to core")
        else:
            logger.info(f"✓ CPU-only configuration (no automatic task mapping)")


def auto_configure_for_partition(partition_name: str) -> HardwareConfig:
    """
    Convenience function to auto-detect hardware for a specific partition.
    
    Args:
        partition_name: Name of the SLURM partition
        
    Returns:
        HardwareConfig with complete hardware info
        
    Example:
        >>> config = auto_configure_for_partition("boost_usr_prod")
        >>> print(config.get_ppr_mapping())
        ppr:4:node
    """
    detector = HardwareDetector()
    return detector.detect_hardware(partition_name=partition_name)


def auto_configure_resource_config(resource_config, partition_name: Optional[str] = None):
    """
    Automatically configure a ResourceConfig with detected hardware.
    
    Args:
        resource_config: ResourceConfig instance to configure
        partition_name: Optional partition name for partition-specific detection
        
    Example:
        >>> from core.config import ResourceConfig
        >>> rc = ResourceConfig()
        >>> auto_configure_resource_config(rc, "boost_usr_prod")
        >>> print(rc.actual_mpi_tasks)  # 4 on Leonardo Booster
        >>> print(rc.cores_per_task)    # 8 on Leonardo Booster
    """
    detector = HardwareDetector()
    detector.detect_hardware(partition_name=partition_name)
    detector.configure_resource_config(resource_config)


if __name__ == '__main__':
    # Test hardware detection
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s - %(message)s'
    )
    
    print("\n" + "="*70)
    print(" Automatic Hardware Detection Test")
    print("="*70 + "\n")
    
    detector = HardwareDetector()
    config = detector.detect_hardware()
    
    print("\n" + "="*70)
    print(config)
    print("="*70)
    
    if config.gpus_per_node > 0:
        print(f"\n✓ System is GPU-accelerated")
        print(f"\nExample mpirun command for 2 nodes:")
        total_tasks = config.mpi_tasks_per_node * 2
        print(f"  mpirun -np {total_tasks} --map-by {config.get_ppr_mapping()} ./your_app")
        print(f"\nThis will run:")
        print(f"  - {config.mpi_tasks_per_node} MPI tasks per node (1 per GPU)")
        print(f"  - {config.cores_per_task} CPU cores per MPI task")
        print(f"  - Total: {total_tasks} MPI tasks across 2 nodes")
    else:
        print(f"\n✓ System is CPU-only (no GPUs detected)")
        print(f"  Would use {config.mpi_tasks_per_node} MPI tasks per node")
