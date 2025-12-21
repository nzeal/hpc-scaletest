#!/usr/bin/env python3
"""
Partition-Aware Hardware Topology Detection

This module provides automatic detection of hardware topology for SLURM partitions
and generates the correct mpirun commands without hardcoded values.

Detection Flow:
1. SLURM environment variables (inside job)
2. sinfo query (login node, pre-submission)
3. scontrol show partition/node (fallback)

Derived Values:
- ranks_per_node = gpus_per_node (for GPU jobs)
- cores_per_rank = cpu_cores_per_node / ranks_per_node

Generated mpirun syntax:
    mpirun -np N --map-by ppr:R:node:PE=C --report-bindings ./bind.sh <executable>

Where:
- N = num_nodes × ranks_per_node
- R = ranks_per_node (GPUs per node)
- C = cores_per_rank

Author: HPC-ScaleTest Contributors
"""

import os
import subprocess
import logging
import re
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Any

logger = logging.getLogger(__name__)


@dataclass
class PartitionTopology:
    """
    Hardware topology for a SLURM partition.
    
    All values are detected at runtime - no hardcoded defaults.
    
    Attributes:
        partition: SLURM partition name
        cpu_cores_per_node: Total CPU cores per node
        gpus_per_node: Number of GPUs per node
        memory_gb: Memory per node in GB
        gpu_model: GPU model name (e.g., "a100")
        detection_method: How topology was detected
        
    Derived (computed automatically):
        ranks_per_node: MPI ranks per node (= GPUs for GPU jobs)
        cores_per_rank: CPU cores per MPI rank
    """
    partition: str
    cpu_cores_per_node: int
    gpus_per_node: int
    memory_gb: float = 0.0
    gpu_model: str = ""
    detection_method: str = ""
    
    # Derived values
    ranks_per_node: int = field(init=False)
    cores_per_rank: int = field(init=False)
    
    def __post_init__(self):
        """Compute derived MPI mapping values."""
        if self.gpus_per_node > 0:
            # GPU job: 1 MPI rank per GPU
            self.ranks_per_node = self.gpus_per_node
            self.cores_per_rank = self.cpu_cores_per_node // self.gpus_per_node
        else:
            # CPU-only job: 1 rank per core
            self.ranks_per_node = self.cpu_cores_per_node
            self.cores_per_rank = 1
    
    def __str__(self) -> str:
        return (
            f"PartitionTopology({self.partition}):\n"
            f"  Detected: {self.cpu_cores_per_node} CPU cores, {self.gpus_per_node} GPUs\n"
            f"  Derived: {self.ranks_per_node} ranks/node, {self.cores_per_rank} cores/rank\n"
            f"  Method: {self.detection_method}"
        )


@dataclass
class JobResources:
    """
    Complete job resource configuration for N nodes.
    
    Contains both SLURM directives (full node allocation) and
    mpirun parameters (actual task distribution).
    """
    num_nodes: int
    topology: PartitionTopology
    
    # SLURM directives (full node resources)
    slurm_ntasks_per_node: int = field(init=False)
    slurm_gres: str = field(init=False)
    
    # MPI parameters (actual distribution)
    mpi_np: int = field(init=False)
    mpi_ranks_per_node: int = field(init=False)
    mpi_cores_per_rank: int = field(init=False)
    
    def __post_init__(self):
        """Compute all resource values."""
        # SLURM: Request full node resources
        self.slurm_ntasks_per_node = self.topology.cpu_cores_per_node
        self.slurm_gres = f"gpu:{self.topology.gpus_per_node}" if self.topology.gpus_per_node > 0 else ""
        
        # MPI: Actual task distribution
        self.mpi_np = self.num_nodes * self.topology.ranks_per_node
        self.mpi_ranks_per_node = self.topology.ranks_per_node
        self.mpi_cores_per_rank = self.topology.cores_per_rank
    
    def get_slurm_directives(self) -> Dict[str, str]:
        """Get SLURM directives as dictionary."""
        directives = {
            'nodes': str(self.num_nodes),
            'partition': self.topology.partition,
            'ntasks-per-node': str(self.slurm_ntasks_per_node),
        }
        if self.slurm_gres:
            directives['gres'] = self.slurm_gres
        return directives
    
    def get_mpirun_command(
        self,
        executable: str,
        args: List[str] = None,
        bind_script: str = "./bind.sh",
        report_bindings: bool = True,
    ) -> str:
        """
        Generate mpirun command string.
        
        Format:
            mpirun -np N --map-by ppr:R:node:PE=C --report-bindings ./bind.sh <exec> <args>
        
        Where:
            N = total MPI ranks = num_nodes × ranks_per_node
            R = ranks_per_node = GPUs per node  
            C = cores_per_rank = CPU cores / ranks
        """
        args = args or []
        
        parts = ['mpirun']
        
        # Total MPI ranks
        parts.append(f'-np {self.mpi_np}')
        
        # Process mapping with PE (Processing Elements) for cores per rank
        parts.append(f'--map-by ppr:{self.mpi_ranks_per_node}:node:PE={self.mpi_cores_per_rank}')
        
        # Report bindings for debugging
        if report_bindings:
            parts.append('--report-bindings')
        
        # GPU binding wrapper
        if self.topology.gpus_per_node > 0 and bind_script:
            parts.append(bind_script)
        
        # Executable and arguments
        parts.append(executable)
        parts.extend(args)
        
        return ' '.join(parts)
    
    def get_mpirun_command_list(
        self,
        executable: str,
        args: List[str] = None,
        bind_script: str = "./bind.sh",
        report_bindings: bool = True,
    ) -> List[str]:
        """Get mpirun command as list of strings."""
        args = args or []
        
        cmd = ['mpirun', '-np', str(self.mpi_np)]
        cmd.extend(['--map-by', f'ppr:{self.mpi_ranks_per_node}:node:PE={self.mpi_cores_per_rank}'])
        
        if report_bindings:
            cmd.append('--report-bindings')
        
        if self.topology.gpus_per_node > 0 and bind_script:
            cmd.append(bind_script)
        
        cmd.append(executable)
        cmd.extend(args)
        
        return cmd


class PartitionDetector:
    """
    Detect hardware topology for SLURM partitions.
    
    Detection Methods (priority order):
    1. SLURM environment variables (if inside a job)
    2. sinfo -p <partition> query
    3. scontrol show partition/node
    
    Usage:
        detector = PartitionDetector()
        topology = detector.detect("boost_usr_prod")
        
        # For Leonardo Booster, this returns:
        # PartitionTopology(cpu_cores_per_node=32, gpus_per_node=4)
        # Derived: ranks_per_node=4, cores_per_rank=8
    """
    
    def __init__(self):
        self._cache: Dict[str, PartitionTopology] = {}
    
    def detect(self, partition: str, force_refresh: bool = False) -> PartitionTopology:
        """
        Detect topology for a SLURM partition.
        
        Args:
            partition: SLURM partition name (e.g., "boost_usr_prod")
            force_refresh: Force re-detection even if cached
        
        Returns:
            PartitionTopology with detected hardware configuration
        
        Raises:
            RuntimeError: If detection fails completely
        """
        if not force_refresh and partition in self._cache:
            return self._cache[partition]
        
        topology = None
        detection_method = ""
        
        # Method 1: SLURM environment (inside job)
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
        logger.info(f"  → Ranks per node: {topology.ranks_per_node}")
        logger.info(f"  → Cores per rank: {topology.cores_per_rank}")
        logger.info(f"  Detection method: {detection_method}")
        
        return topology
    
    def _is_inside_slurm_job(self, partition: str) -> bool:
        """Check if running inside a SLURM job for this partition."""
        return (
            os.environ.get('SLURM_JOB_ID') is not None and
            os.environ.get('SLURM_JOB_PARTITION') == partition
        )
    
    def _detect_from_slurm_env(self, partition: str) -> Optional[PartitionTopology]:
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
            
            return PartitionTopology(
                partition=partition,
                cpu_cores_per_node=cpu_cores,
                gpus_per_node=gpus,
            )
        except (ValueError, KeyError) as e:
            logger.debug(f"SLURM env detection failed: {e}")
            return None
    
    def _detect_from_sinfo(self, partition: str) -> Optional[PartitionTopology]:
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
            
            # Parse first line
            line = result.stdout.strip().split('\n')[0]
            parts = line.split()
            
            cpu_cores = None
            gpus = 0
            memory_gb = 0.0
            gpu_model = ""
            
            for part in parts:
                # CPU cores (first numeric)
                if cpu_cores is None and part.isdigit():
                    cpu_cores = int(part)
                
                # GPU info: gpu:type:count or gpu:count or (null)
                elif 'gpu' in part.lower() and part != '(null)':
                    match = re.search(r'gpu:(?:([^:]+):)?(\d+)', part, re.IGNORECASE)
                    if match:
                        gpu_model = match.group(1) or ""
                        gpus = int(match.group(2))
                
                # Memory (second numeric, in MB)
                elif part.isdigit() and cpu_cores is not None:
                    memory_gb = int(part) / 1024.0
            
            if cpu_cores is None:
                return None
            
            return PartitionTopology(
                partition=partition,
                cpu_cores_per_node=cpu_cores,
                gpus_per_node=gpus,
                memory_gb=memory_gb,
                gpu_model=gpu_model,
            )
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.debug(f"sinfo detection failed: {e}")
            return None
    
    def _detect_from_scontrol(self, partition: str) -> Optional[PartitionTopology]:
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
            
            # If no CPU info, query a node
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
            
            return PartitionTopology(
                partition=partition,
                cpu_cores_per_node=cpu_cores,
                gpus_per_node=gpus,
            )
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.debug(f"scontrol detection failed: {e}")
            return None


# =============================================================================
# Module-level singleton
# =============================================================================

_detector: Optional[PartitionDetector] = None

def get_partition_detector() -> PartitionDetector:
    """Get or create the global partition detector."""
    global _detector
    if _detector is None:
        _detector = PartitionDetector()
    return _detector


def detect_topology(partition: str) -> PartitionTopology:
    """
    Convenience function to detect partition topology.
    
    Args:
        partition: SLURM partition name
    
    Returns:
        PartitionTopology with detected configuration
    
    Example:
        >>> topo = detect_topology("boost_usr_prod")
        >>> print(f"{topo.cpu_cores_per_node} cores, {topo.gpus_per_node} GPUs")
        32 cores, 4 GPUs
        >>> print(f"{topo.ranks_per_node} ranks/node, {topo.cores_per_rank} cores/rank")
        4 ranks/node, 8 cores/rank
    """
    return get_partition_detector().detect(partition)


def generate_job_resources(partition: str, num_nodes: int) -> JobResources:
    """
    Generate complete job resource configuration.
    
    Args:
        partition: SLURM partition name
        num_nodes: Number of nodes
    
    Returns:
        JobResources with SLURM directives and mpirun parameters
    
    Example:
        >>> res = generate_job_resources("boost_usr_prod", 4)
        >>> print(res.get_mpirun_command("$BINARY/iPIC3D", ["os-stdin"]))
        mpirun -np 16 --map-by ppr:4:node:PE=8 --report-bindings ./bind.sh $BINARY/iPIC3D os-stdin
    """
    topology = detect_topology(partition)
    return JobResources(num_nodes=num_nodes, topology=topology)


def generate_bind_script(output_path: str = "./bind.sh") -> str:
    """
    Generate GPU binding script using OMPI_COMM_WORLD_LOCAL_RANK.
    
    Args:
        output_path: Path to write script
    
    Returns:
        Path to generated script
    """
    script = '''#!/bin/bash
# GPU Binding Script - Generated by HPC-ScaleTest
# Binds each MPI rank to a unique GPU via CUDA_VISIBLE_DEVICES

# Determine local rank from environment
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
    
    with open(output_path, 'w') as f:
        f.write(script)
    os.chmod(output_path, 0o755)
    
    return output_path


# =============================================================================
# Self-test
# =============================================================================

if __name__ == '__main__':
    print("=" * 70)
    print(" Partition Topology Detection - Self Test")
    print("=" * 70)
    
    # Test with mock data (simulating Leonardo Booster)
    print("\n[Test 1] Leonardo Booster: 32 cores, 4 GPUs")
    
    topo = PartitionTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        detection_method="mock"
    )
    
    print(f"  Detected: {topo.cpu_cores_per_node} CPU cores, {topo.gpus_per_node} GPUs")
    print(f"  Derived:  {topo.ranks_per_node} ranks/node, {topo.cores_per_rank} cores/rank")
    
    # Test mpirun commands for different node counts
    for N in [1, 2, 4]:
        res = JobResources(num_nodes=N, topology=topo)
        cmd = res.get_mpirun_command("$BINARY/iPIC3D", ["os-stdin"])
        
        print(f"\n  {N} node(s):")
        print(f"    SLURM: --ntasks-per-node={res.slurm_ntasks_per_node} --gres={res.slurm_gres}")
        print(f"    {cmd}")
        
        # Verify formula: -np = N × GPUs
        expected_np = N * topo.gpus_per_node
        assert f'-np {expected_np}' in cmd, f"Expected -np {expected_np}"
        assert f'ppr:{topo.gpus_per_node}:node:PE={topo.cores_per_rank}' in cmd
    
    # Test LUMI: 128 cores, 8 GPUs
    print("\n" + "-" * 70)
    print("\n[Test 2] LUMI: 128 cores, 8 GPUs")
    
    lumi = PartitionTopology(
        partition="gpu",
        cpu_cores_per_node=128,
        gpus_per_node=8,
        detection_method="mock"
    )
    
    print(f"  Detected: {lumi.cpu_cores_per_node} CPU cores, {lumi.gpus_per_node} GPUs")
    print(f"  Derived:  {lumi.ranks_per_node} ranks/node, {lumi.cores_per_rank} cores/rank")
    
    res = JobResources(num_nodes=4, topology=lumi)
    cmd = res.get_mpirun_command("./app", ["input.dat"])
    
    print(f"\n  4 nodes:")
    print(f"    {cmd}")
    
    assert '-np 32' in cmd  # 4 nodes × 8 GPUs
    assert 'ppr:8:node:PE=16' in cmd  # 8 ranks, 128/8=16 cores each
    
    print("\n" + "=" * 70)
    print(" All tests PASSED")
    print("=" * 70)
