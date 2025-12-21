#!/usr/bin/env python3
"""
core/hardware.py - Unified Hardware Topology Detection

This is the SINGLE SOURCE OF TRUTH for hardware topology detection in HPC-ScaleTest.
All other modules should import from here.

Design Principles:
1. NO HARDCODED VALUES - All topology detected at runtime
2. PARTITION-AWARE - Queries SLURM for partition-specific hardware
3. SINGLE RESPONSIBILITY - Only topology detection, no job script generation

Detection Priority:
1. SLURM environment variables (inside job)
2. sinfo/scontrol queries (login node)
3. System introspection (/proc, nvidia-smi)

Key Concepts:
- cpu_cores_per_node: Total CPU cores available on a node
- gpus_per_node: Number of GPUs per node
- mpi_ranks_per_node: Number of MPI processes per node (= GPUs for GPU jobs)
- cores_per_rank: CPU cores allocated to each MPI rank

For GPU jobs:
    mpi_ranks_per_node = gpus_per_node
    cores_per_rank = cpu_cores_per_node / gpus_per_node

For CPU jobs:
    mpi_ranks_per_node = cpu_cores_per_node
    cores_per_rank = 1

Author: HPC-ScaleTest Contributors
"""

import os
import subprocess
import logging
import re
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple
from enum import Enum

logger = logging.getLogger(__name__)


class GPUVendor(Enum):
    """GPU vendor types."""
    NVIDIA = "nvidia"
    AMD = "amd"
    INTEL = "intel"
    NONE = "none"


@dataclass
class HardwareTopology:
    """
    Hardware topology for a compute node.
    
    This is the canonical representation of node hardware configuration.
    All values are detected at runtime - no hardcoded defaults.
    
    Attributes:
        partition: SLURM partition name
        cpu_cores_per_node: Total CPU cores per node
        gpus_per_node: Number of GPUs per node
        memory_gb: Memory per node in GB
        gpu_vendor: GPU vendor (nvidia, amd, intel, none)
        gpu_model: GPU model name (e.g., "a100")
        detection_method: How topology was detected
        
    Derived (computed in __post_init__):
        mpi_ranks_per_node: MPI processes per node
        cores_per_rank: CPU cores per MPI process
    """
    partition: str
    cpu_cores_per_node: int
    gpus_per_node: int = 0
    memory_gb: float = 0.0
    gpu_vendor: GPUVendor = GPUVendor.NONE
    gpu_model: str = ""
    detection_method: str = ""
    
    # Derived values (computed automatically)
    mpi_ranks_per_node: int = field(init=False)
    cores_per_rank: int = field(init=False)
    
    def __post_init__(self):
        """Compute derived MPI mapping values."""
        if self.gpus_per_node > 0:
            # GPU job: 1 MPI rank per GPU
            self.mpi_ranks_per_node = self.gpus_per_node
            self.cores_per_rank = self.cpu_cores_per_node // self.gpus_per_node
        else:
            # CPU-only job: 1 MPI rank per core
            self.mpi_ranks_per_node = self.cpu_cores_per_node
            self.cores_per_rank = 1
    
    def is_gpu_partition(self) -> bool:
        """Check if this is a GPU partition."""
        return self.gpus_per_node > 0
    
    def total_mpi_ranks(self, num_nodes: int) -> int:
        """Calculate total MPI ranks for N nodes."""
        return num_nodes * self.mpi_ranks_per_node
    
    def __str__(self) -> str:
        return (
            f"HardwareTopology({self.partition}):\n"
            f"  Detected: {self.cpu_cores_per_node} CPU cores, {self.gpus_per_node} GPUs\n"
            f"  Derived: {self.mpi_ranks_per_node} MPI ranks/node, {self.cores_per_rank} cores/rank\n"
            f"  Method: {self.detection_method}"
        )


@dataclass
class MPIConfiguration:
    """
    MPI launch configuration derived from topology and node count.
    
    This class encapsulates all values needed for MPI execution on N nodes.
    """
    num_nodes: int
    topology: HardwareTopology
    
    # Total values for the job
    total_mpi_ranks: int = field(init=False)
    mpi_ranks_per_node: int = field(init=False)
    cores_per_rank: int = field(init=False)
    total_gpus: int = field(init=False)
    
    def __post_init__(self):
        """Compute configuration values."""
        self.mpi_ranks_per_node = self.topology.mpi_ranks_per_node
        self.cores_per_rank = self.topology.cores_per_rank
        self.total_mpi_ranks = self.num_nodes * self.mpi_ranks_per_node
        self.total_gpus = self.num_nodes * self.topology.gpus_per_node
    
    def get_mpirun_args(self, report_bindings: bool = True) -> List[str]:
        """
        Generate mpirun arguments for OpenMPI.
        
        Returns arguments for: mpirun -np N --map-by ppr:R:node:PE=C
        
        Args:
            report_bindings: Include --report-bindings flag
            
        Returns:
            List of mpirun arguments (not including mpirun itself)
        """
        args = [
            '-np', str(self.total_mpi_ranks),
            '--map-by', f'ppr:{self.mpi_ranks_per_node}:node:PE={self.cores_per_rank}',
        ]
        
        if report_bindings:
            args.append('--report-bindings')
        
        return args
    
    def get_mpirun_command_string(self, report_bindings: bool = True) -> str:
        """Get mpirun arguments as a string."""
        args = self.get_mpirun_args(report_bindings)
        return ' '.join(args)


class TopologyDetector:
    """
    Detect hardware topology for SLURM partitions.
    
    Detection Methods (priority order):
    1. SLURM environment variables (if inside a job)
    2. sinfo query for partition hardware
    3. scontrol query for partition/node info
    4. System introspection (fallback)
    
    Usage:
        detector = TopologyDetector()
        topology = detector.detect("boost_usr_prod")
        
        # For Leonardo Booster, returns:
        # HardwareTopology(cpu_cores=32, gpus=4)
        # Derived: mpi_ranks=4, cores_per_rank=8
    """
    
    def __init__(self):
        self._cache: Dict[str, HardwareTopology] = {}
    
    def detect(self, partition: str, force_refresh: bool = False) -> HardwareTopology:
        """
        Detect hardware topology for a SLURM partition.
        
        Args:
            partition: SLURM partition name (e.g., "boost_usr_prod")
            force_refresh: Force re-detection even if cached
            
        Returns:
            HardwareTopology with detected configuration
            
        Raises:
            RuntimeError: If detection fails completely
        """
        if not partition:
            raise ValueError("Partition name is required for topology detection")
        
        if not force_refresh and partition in self._cache:
            logger.debug(f"Using cached topology for partition {partition}")
            return self._cache[partition]
        
        topology = None
        detection_method = ""
        
        # Method 1: SLURM environment variables (inside job)
        if self._is_inside_slurm_job(partition):
            topology = self._detect_from_slurm_env(partition)
            if topology:
                detection_method = "SLURM environment"
        
        # Method 2: sinfo query
        if topology is None:
            topology = self._detect_from_sinfo(partition)
            if topology:
                detection_method = "sinfo query"
        
        # Method 3: scontrol query
        if topology is None:
            topology = self._detect_from_scontrol(partition)
            if topology:
                detection_method = "scontrol query"
        
        # Method 4: System introspection (fallback)
        if topology is None:
            topology = self._detect_from_system(partition)
            if topology:
                detection_method = "system introspection"
        
        if topology is None:
            raise RuntimeError(
                f"Failed to detect topology for partition '{partition}'. "
                f"Ensure SLURM is available and partition exists."
            )
        
        topology.detection_method = detection_method
        self._cache[partition] = topology
        
        logger.info(f"Detected topology for '{partition}':")
        logger.info(f"  CPU cores per node: {topology.cpu_cores_per_node}")
        logger.info(f"  GPUs per node: {topology.gpus_per_node}")
        logger.info(f"  MPI ranks per node: {topology.mpi_ranks_per_node}")
        logger.info(f"  Cores per rank: {topology.cores_per_rank}")
        logger.info(f"  Detection method: {detection_method}")
        
        return topology
    
    def get_mpi_config(self, partition: str, num_nodes: int) -> MPIConfiguration:
        """
        Get complete MPI configuration for a job.
        
        Args:
            partition: SLURM partition name
            num_nodes: Number of nodes for the job
            
        Returns:
            MPIConfiguration with all computed values
        """
        topology = self.detect(partition)
        return MPIConfiguration(num_nodes=num_nodes, topology=topology)
    
    def _is_inside_slurm_job(self, partition: str) -> bool:
        """Check if running inside a SLURM job for this partition."""
        return (
            os.environ.get('SLURM_JOB_ID') is not None and
            os.environ.get('SLURM_JOB_PARTITION') == partition
        )
    
    def _detect_from_slurm_env(self, partition: str) -> Optional[HardwareTopology]:
        """Detect from SLURM environment variables."""
        try:
            cpu_cores = None
            gpus = 0
            
            # CPU cores
            if 'SLURM_CPUS_ON_NODE' in os.environ:
                cpu_cores = int(os.environ['SLURM_CPUS_ON_NODE'])
            elif 'SLURM_JOB_CPUS_PER_NODE' in os.environ:
                # Format: "32" or "32(x4)"
                val = os.environ['SLURM_JOB_CPUS_PER_NODE']
                cpu_cores = int(val.split('(')[0])
            
            if cpu_cores is None:
                return None
            
            # GPUs
            if 'SLURM_GPUS_ON_NODE' in os.environ:
                gpus = int(os.environ['SLURM_GPUS_ON_NODE'])
            elif 'SLURM_GPUS_PER_NODE' in os.environ:
                gpus = int(os.environ['SLURM_GPUS_PER_NODE'])
            elif 'CUDA_VISIBLE_DEVICES' in os.environ:
                devs = os.environ['CUDA_VISIBLE_DEVICES']
                if devs and devs.lower() not in ('', 'nodevfiles'):
                    gpus = len([d for d in devs.split(',') if d.strip()])
            
            return HardwareTopology(
                partition=partition,
                cpu_cores_per_node=cpu_cores,
                gpus_per_node=gpus,
            )
        except (ValueError, KeyError) as e:
            logger.debug(f"SLURM env detection failed: {e}")
            return None
    
    def _detect_from_sinfo(self, partition: str) -> Optional[HardwareTopology]:
        """Detect from sinfo query."""
        try:
            # Query: CPUs, GRES (GPUs), Memory
            result = subprocess.run(
                ['sinfo', '-p', partition, '-N', '-h', '-o', '%c %G %m'],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0 or not result.stdout.strip():
                return None
            
            # Parse first line (representative node)
            line = result.stdout.strip().split('\n')[0]
            parts = line.split()
            
            cpu_cores = None
            gpus = 0
            memory_gb = 0.0
            gpu_model = ""
            gpu_vendor = GPUVendor.NONE
            
            for part in parts:
                # CPU cores (first numeric value)
                if cpu_cores is None and part.isdigit():
                    cpu_cores = int(part)
                
                # GPU info: gpu:type:count or gpu:count or (null)
                elif 'gpu' in part.lower() and part != '(null)':
                    match = re.search(r'gpu:(?:([^:]+):)?(\d+)', part, re.IGNORECASE)
                    if match:
                        gpu_model = match.group(1) or ""
                        gpus = int(match.group(2))
                        
                        # Detect vendor from model
                        if any(x in gpu_model.lower() for x in ['a100', 'v100', 'h100', 'a30', 'a40']):
                            gpu_vendor = GPUVendor.NVIDIA
                        elif any(x in gpu_model.lower() for x in ['mi100', 'mi200', 'mi250', 'mi300']):
                            gpu_vendor = GPUVendor.AMD
                
                # Memory (second numeric, in MB)
                elif part.isdigit() and cpu_cores is not None:
                    memory_gb = int(part) / 1024.0
            
            if cpu_cores is None:
                return None
            
            return HardwareTopology(
                partition=partition,
                cpu_cores_per_node=cpu_cores,
                gpus_per_node=gpus,
                memory_gb=memory_gb,
                gpu_vendor=gpu_vendor,
                gpu_model=gpu_model,
            )
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.debug(f"sinfo detection failed: {e}")
            return None
    
    def _detect_from_scontrol(self, partition: str) -> Optional[HardwareTopology]:
        """Detect from scontrol show partition/node."""
        try:
            # Query partition
            result = subprocess.run(
                ['scontrol', 'show', 'partition', partition],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0:
                return None
            
            output = result.stdout
            cpu_cores = None
            gpus = 0
            
            # MaxCPUsPerNode
            match = re.search(r'MaxCPUsPerNode=(\d+)', output)
            if match:
                cpu_cores = int(match.group(1))
            
            # GRES for GPUs
            match = re.search(r'GRES=([^\s]+)', output)
            if match:
                gres = match.group(1)
                gpu_match = re.search(r'gpu:(?:\w+:)?(\d+)', gres)
                if gpu_match:
                    gpus = int(gpu_match.group(1))
            
            # If no CPU info from partition, query a node
            if cpu_cores is None:
                match = re.search(r'Nodes=([^\s,\[]+)', output)
                if match:
                    node = match.group(1)
                    node_result = subprocess.run(
                        ['scontrol', 'show', 'node', node],
                        capture_output=True,
                        text=True,
                        timeout=10
                    )
                    if node_result.returncode == 0:
                        cpu_match = re.search(r'CPUTot=(\d+)', node_result.stdout)
                        if cpu_match:
                            cpu_cores = int(cpu_match.group(1))
                        
                        if gpus == 0:
                            gres_match = re.search(r'Gres=([^\s]+)', node_result.stdout)
                            if gres_match:
                                gres = gres_match.group(1)
                                gpu_match = re.search(r'gpu:(?:\w+:)?(\d+)', gres)
                                if gpu_match:
                                    gpus = int(gpu_match.group(1))
            
            if cpu_cores is None:
                return None
            
            return HardwareTopology(
                partition=partition,
                cpu_cores_per_node=cpu_cores,
                gpus_per_node=gpus,
            )
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.debug(f"scontrol detection failed: {e}")
            return None
    
    def _detect_from_system(self, partition: str) -> Optional[HardwareTopology]:
        """Detect from system introspection (fallback)."""
        try:
            # Get CPU count
            cpu_cores = os.cpu_count() or 1
            
            # Try to get physical cores from lscpu
            try:
                result = subprocess.run(
                    ['lscpu'],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0:
                    # Look for "Core(s) per socket" and "Socket(s)"
                    cores_match = re.search(r'Core\(s\) per socket:\s+(\d+)', result.stdout)
                    sockets_match = re.search(r'Socket\(s\):\s+(\d+)', result.stdout)
                    if cores_match and sockets_match:
                        cpu_cores = int(cores_match.group(1)) * int(sockets_match.group(1))
            except (subprocess.TimeoutExpired, FileNotFoundError):
                pass
            
            # Detect GPUs via nvidia-smi
            gpus = 0
            gpu_vendor = GPUVendor.NONE
            gpu_model = ""
            
            try:
                result = subprocess.run(
                    ['nvidia-smi', '--query-gpu=count', '--format=csv,noheader'],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0:
                    gpus = int(result.stdout.strip().split('\n')[0])
                    gpu_vendor = GPUVendor.NVIDIA
                    
                    # Get model
                    model_result = subprocess.run(
                        ['nvidia-smi', '--query-gpu=name', '--format=csv,noheader'],
                        capture_output=True,
                        text=True,
                        timeout=5
                    )
                    if model_result.returncode == 0:
                        gpu_model = model_result.stdout.strip().split('\n')[0]
            except (subprocess.TimeoutExpired, FileNotFoundError):
                pass
            
            # Try AMD GPUs
            if gpus == 0:
                try:
                    result = subprocess.run(
                        ['rocm-smi', '--showproductname'],
                        capture_output=True,
                        text=True,
                        timeout=5
                    )
                    if result.returncode == 0:
                        # Count GPU entries
                        gpus = result.stdout.count('GPU')
                        gpu_vendor = GPUVendor.AMD
                except (subprocess.TimeoutExpired, FileNotFoundError):
                    pass
            
            return HardwareTopology(
                partition=partition,
                cpu_cores_per_node=cpu_cores,
                gpus_per_node=gpus,
                gpu_vendor=gpu_vendor,
                gpu_model=gpu_model,
            )
        except Exception as e:
            logger.debug(f"System introspection failed: {e}")
            return None


# =============================================================================
# Module-level singleton and convenience functions
# =============================================================================

_detector: Optional[TopologyDetector] = None


def get_topology_detector() -> TopologyDetector:
    """Get or create the global topology detector."""
    global _detector
    if _detector is None:
        _detector = TopologyDetector()
    return _detector


def detect_topology(partition: str) -> HardwareTopology:
    """
    Convenience function to detect partition topology.
    
    Args:
        partition: SLURM partition name
        
    Returns:
        HardwareTopology with detected configuration
        
    Example:
        >>> topo = detect_topology("boost_usr_prod")
        >>> print(f"{topo.cpu_cores_per_node} cores, {topo.gpus_per_node} GPUs")
        32 cores, 4 GPUs
    """
    return get_topology_detector().detect(partition)


def get_mpi_config(partition: str, num_nodes: int) -> MPIConfiguration:
    """
    Convenience function to get MPI configuration.
    
    Args:
        partition: SLURM partition name
        num_nodes: Number of nodes
        
    Returns:
        MPIConfiguration with all computed values
        
    Example:
        >>> cfg = get_mpi_config("boost_usr_prod", 4)
        >>> print(cfg.get_mpirun_args())
        ['-np', '16', '--map-by', 'ppr:4:node:PE=8', '--report-bindings']
    """
    return get_topology_detector().get_mpi_config(partition, num_nodes)


# =============================================================================
# GPU Binding Script Generation
# =============================================================================

GPU_BIND_SCRIPT = '''#!/bin/bash
# GPU Binding Script - Generated by HPC-ScaleTest
# Binds each MPI rank to a unique GPU via CUDA_VISIBLE_DEVICES

# Determine local rank from available environment variables
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

export CUDA_VISIBLE_DEVICES=$LOCAL_RANK
exec "$@"
'''


def generate_gpu_bind_script(output_path: str = "./bind.sh") -> str:
    """
    Generate GPU binding script.
    
    Args:
        output_path: Path to write the script
        
    Returns:
        Path to generated script
    """
    with open(output_path, 'w') as f:
        f.write(GPU_BIND_SCRIPT)
    os.chmod(output_path, 0o755)
    logger.info(f"Generated GPU binding script: {output_path}")
    return output_path


# =============================================================================
# Self-test
# =============================================================================

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
    
    print("=" * 70)
    print(" Hardware Topology Detection - Self Test")
    print("=" * 70)
    
    # Test with mock topology (simulating Leonardo Booster)
    print("\n[Test] Leonardo Booster: 32 cores, 4 GPUs")
    
    topo = HardwareTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        gpu_vendor=GPUVendor.NVIDIA,
        gpu_model="a100",
        detection_method="test"
    )
    
    print(f"  {topo}")
    assert topo.mpi_ranks_per_node == 4, "Expected 4 ranks/node"
    assert topo.cores_per_rank == 8, "Expected 8 cores/rank"
    print("  ✓ Derived values correct")
    
    # Test MPI configuration
    for N in [1, 2, 4]:
        cfg = MPIConfiguration(num_nodes=N, topology=topo)
        args = cfg.get_mpirun_args()
        cmd_str = ' '.join(['mpirun'] + args)
        
        print(f"\n  {N} node(s):")
        print(f"    {cmd_str}")
        
        expected_np = N * 4
        assert f'-np {expected_np}' in cmd_str
        assert 'ppr:4:node:PE=8' in cmd_str
    
    print("\n" + "=" * 70)
    print(" All tests PASSED")
    print("=" * 70)
