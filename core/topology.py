#!/usr/bin/env python3
"""
Centralized Hardware Topology Detection Module

This module provides the SINGLE SOURCE OF TRUTH for hardware topology detection
in HPC-ScaleTest. All other modules should import from here rather than
implementing their own detection logic.

Design Principles:
1. NO HARDCODED VALUES - All topology is detected at runtime
2. SLURM-FIRST - Use SLURM environment variables when available
3. FALLBACK CHAIN - Multiple detection methods with graceful degradation
4. EXPLICIT ASSUMPTIONS - All assumptions are documented and configurable

Detection Strategy:
1. Job Runtime (compute node): SLURM env vars → system introspection
2. Job Submission (login node): sinfo/scontrol → cached partition data

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


class DetectionContext(Enum):
    """Context in which topology detection is running."""
    COMPUTE_NODE = "compute"  # Inside a SLURM job on compute node
    LOGIN_NODE = "login"      # On login node, querying for job preparation
    UNKNOWN = "unknown"       # Cannot determine context


class GPUVendor(Enum):
    """Supported GPU vendors."""
    NVIDIA = "nvidia"
    AMD = "amd"
    INTEL = "intel"
    NONE = "none"


@dataclass
class NodeTopology:
    """
    Complete hardware topology for a single node.
    
    This is the canonical representation of node hardware.
    All values are detected at runtime - no defaults.
    
    Attributes:
        cpu_cores: Total CPU cores available on the node
        cpu_sockets: Number of CPU sockets
        cores_per_socket: Cores per CPU socket
        threads_per_core: Hardware threads per core (hyperthreading)
        physical_cores: Physical cores (cpu_cores / threads_per_core)
        gpus: Number of GPUs on the node
        gpu_vendor: GPU vendor (nvidia, amd, intel, none)
        gpu_model: GPU model name if detected
        gpu_memory_gb: GPU memory in GB if detected
        memory_gb: Total system memory in GB
        partition: SLURM partition name (if known)
        detection_method: How the topology was detected
    """
    cpu_cores: int
    cpu_sockets: int = 1
    cores_per_socket: int = 0
    threads_per_core: int = 1
    physical_cores: int = 0
    gpus: int = 0
    gpu_vendor: GPUVendor = GPUVendor.NONE
    gpu_model: str = ""
    gpu_memory_gb: float = 0.0
    memory_gb: float = 0.0
    partition: str = ""
    detection_method: str = ""
    
    def __post_init__(self):
        """Compute derived values."""
        if self.cores_per_socket == 0 and self.cpu_sockets > 0:
            self.cores_per_socket = self.cpu_cores // self.cpu_sockets
        if self.physical_cores == 0:
            self.physical_cores = self.cpu_cores // self.threads_per_core


@dataclass
class MPIMapping:
    """
    MPI process mapping configuration derived from topology.
    
    This class encapsulates the logic for mapping MPI ranks to hardware
    resources based on the detected topology.
    
    Attributes:
        ranks_per_node: Number of MPI ranks per node
        cores_per_rank: CPU cores allocated to each rank
        gpus_per_rank: GPUs allocated to each rank (typically 0 or 1)
        total_ranks: Total MPI ranks across all nodes (if known)
        binding_strategy: How ranks are bound to resources
    """
    ranks_per_node: int
    cores_per_rank: int
    gpus_per_rank: int = 0
    total_ranks: int = 0
    binding_strategy: str = "core"
    
    def validate(self, topology: NodeTopology) -> Tuple[bool, str]:
        """
        Validate that this mapping is consistent with the topology.
        
        Returns:
            Tuple of (is_valid, error_message)
        """
        # Check CPU allocation
        total_cores_used = self.ranks_per_node * self.cores_per_rank
        if total_cores_used > topology.cpu_cores:
            return False, (
                f"CPU oversubscription: {self.ranks_per_node} ranks × "
                f"{self.cores_per_rank} cores = {total_cores_used} > "
                f"{topology.cpu_cores} available"
            )
        
        # Check GPU allocation
        if self.gpus_per_rank > 0:
            total_gpus_used = self.ranks_per_node * self.gpus_per_rank
            if total_gpus_used > topology.gpus:
                return False, (
                    f"GPU oversubscription: {self.ranks_per_node} ranks × "
                    f"{self.gpus_per_rank} GPUs = {total_gpus_used} > "
                    f"{topology.gpus} available"
                )
        
        # Warn about unused resources
        unused_cores = topology.cpu_cores - total_cores_used
        if unused_cores > 0:
            logger.debug(f"Note: {unused_cores} CPU cores will be unused per node")
        
        return True, ""


class TopologyDetector:
    """
    Unified hardware topology detector.
    
    This class implements a detection chain that tries multiple methods
    to determine hardware topology, with no hardcoded fallback values.
    
    Detection Methods (in order of preference):
    1. SLURM environment variables (most reliable inside a job)
    2. SLURM partition query via sinfo/scontrol (on login node)
    3. Direct system introspection (/proc, nvidia-smi, etc.)
    
    Usage:
        detector = TopologyDetector()
        topology = detector.detect()
        mapping = detector.compute_mpi_mapping(topology, num_nodes=4)
    """
    
    # SLURM environment variables we use for detection
    SLURM_ENV_VARS = {
        'cpus': [
            'SLURM_CPUS_ON_NODE',      # CPUs on this node
            'SLURM_JOB_CPUS_PER_NODE', # CPUs per node in job
            'SLURM_CPUS_PER_TASK',     # CPUs per task
        ],
        'gpus': [
            'SLURM_GPUS_ON_NODE',      # GPUs on this node
            'SLURM_GPUS_PER_NODE',     # GPUs per node requested
            'CUDA_VISIBLE_DEVICES',     # NVIDIA GPU visibility
            'ROCR_VISIBLE_DEVICES',     # AMD GPU visibility
            'ZE_AFFINITY_MASK',         # Intel GPU visibility
        ],
        'job': [
            'SLURM_JOB_ID',            # Indicates we're in a job
            'SLURM_JOB_PARTITION',     # Partition name
            'SLURM_JOB_NUM_NODES',     # Number of nodes
            'SLURM_NODELIST',          # Node list
            'SLURM_TASKS_PER_NODE',    # Tasks per node
        ]
    }
    
    def __init__(self):
        """Initialize the topology detector."""
        self._context: Optional[DetectionContext] = None
        self._cache: Dict[str, NodeTopology] = {}
    
    @property
    def context(self) -> DetectionContext:
        """Determine the current detection context."""
        if self._context is None:
            if os.environ.get('SLURM_JOB_ID'):
                self._context = DetectionContext.COMPUTE_NODE
            elif self._check_slurm_available():
                self._context = DetectionContext.LOGIN_NODE
            else:
                self._context = DetectionContext.UNKNOWN
        return self._context
    
    def _check_slurm_available(self) -> bool:
        """Check if SLURM commands are available."""
        try:
            result = subprocess.run(
                ['sinfo', '--version'],
                capture_output=True,
                timeout=5
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            return False
    
    def detect(self, partition: Optional[str] = None) -> NodeTopology:
        """
        Detect hardware topology.
        
        Args:
            partition: Optional SLURM partition name. If not provided,
                      uses SLURM_JOB_PARTITION or attempts auto-detection.
        
        Returns:
            NodeTopology with detected hardware configuration.
        
        Raises:
            RuntimeError: If topology cannot be detected.
        """
        # Check cache
        cache_key = partition or "_default"
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # Determine partition
        if partition is None:
            partition = os.environ.get('SLURM_JOB_PARTITION', '')
        
        # Try detection methods in order
        topology = None
        methods_tried = []
        
        if self.context == DetectionContext.COMPUTE_NODE:
            # Inside a job - use SLURM env vars + system introspection
            topology = self._detect_from_slurm_env()
            if topology:
                methods_tried.append("SLURM environment")
        
        if topology is None and partition:
            # Try partition query
            topology = self._detect_from_partition(partition)
            if topology:
                methods_tried.append(f"partition query ({partition})")
        
        if topology is None:
            # Fall back to direct system introspection
            topology = self._detect_from_system()
            if topology:
                methods_tried.append("system introspection")
        
        if topology is None:
            raise RuntimeError(
                f"Failed to detect hardware topology. "
                f"Methods tried: {', '.join(methods_tried) or 'none'}. "
                f"Context: {self.context.value}. "
                f"Ensure SLURM is available or run on a compute node."
            )
        
        topology.partition = partition
        topology.detection_method = " + ".join(methods_tried)
        
        # Cache the result
        self._cache[cache_key] = topology
        
        logger.info(f"Detected topology via {topology.detection_method}:")
        logger.info(f"  CPU cores: {topology.cpu_cores}")
        logger.info(f"  GPUs: {topology.gpus} ({topology.gpu_vendor.value})")
        if topology.gpu_model:
            logger.info(f"  GPU model: {topology.gpu_model}")
        
        return topology
    
    def _detect_from_slurm_env(self) -> Optional[NodeTopology]:
        """
        Detect topology from SLURM environment variables.
        
        This is the most reliable method when running inside a SLURM job.
        """
        # Get CPU count
        cpu_cores = None
        for var in self.SLURM_ENV_VARS['cpus']:
            value = os.environ.get(var)
            if value:
                try:
                    # Handle formats like "32" or "32(x2)" 
                    cpu_cores = int(value.split('(')[0])
                    logger.debug(f"CPU cores from {var}: {cpu_cores}")
                    break
                except ValueError:
                    continue
        
        if cpu_cores is None:
            # Try system detection for CPUs
            cpu_cores = self._detect_cpu_count()
        
        if cpu_cores is None:
            return None
        
        # Get GPU count
        gpus = 0
        gpu_vendor = GPUVendor.NONE
        gpu_model = ""
        
        # Check SLURM GPU variables
        slurm_gpus = os.environ.get('SLURM_GPUS_ON_NODE')
        if slurm_gpus:
            try:
                gpus = int(slurm_gpus)
            except ValueError:
                pass
        
        # Check CUDA_VISIBLE_DEVICES
        if gpus == 0:
            cuda_devices = os.environ.get('CUDA_VISIBLE_DEVICES', '')
            if cuda_devices and cuda_devices.lower() not in ('', 'nodevfiles'):
                gpus = len([d for d in cuda_devices.split(',') if d.strip()])
                gpu_vendor = GPUVendor.NVIDIA
        
        # Check ROCR_VISIBLE_DEVICES (AMD)
        if gpus == 0:
            rocr_devices = os.environ.get('ROCR_VISIBLE_DEVICES', '')
            if rocr_devices:
                gpus = len([d for d in rocr_devices.split(',') if d.strip()])
                gpu_vendor = GPUVendor.AMD
        
        # Try nvidia-smi for GPU detection and model
        if gpus == 0 or gpu_vendor == GPUVendor.NONE:
            nvidia_info = self._detect_nvidia_gpus()
            if nvidia_info:
                gpus = nvidia_info['count']
                gpu_vendor = GPUVendor.NVIDIA
                gpu_model = nvidia_info.get('model', '')
        
        # Try rocm-smi for AMD GPUs
        if gpus == 0:
            amd_info = self._detect_amd_gpus()
            if amd_info:
                gpus = amd_info['count']
                gpu_vendor = GPUVendor.AMD
                gpu_model = amd_info.get('model', '')
        
        # Get memory info
        memory_gb = self._detect_memory()
        
        # Detect CPU topology
        cpu_info = self._detect_cpu_topology()
        
        return NodeTopology(
            cpu_cores=cpu_cores,
            cpu_sockets=cpu_info.get('sockets', 1),
            cores_per_socket=cpu_info.get('cores_per_socket', cpu_cores),
            threads_per_core=cpu_info.get('threads_per_core', 1),
            gpus=gpus,
            gpu_vendor=gpu_vendor,
            gpu_model=gpu_model,
            memory_gb=memory_gb,
        )
    
    def _detect_from_partition(self, partition: str) -> Optional[NodeTopology]:
        """
        Detect topology by querying SLURM partition information.
        
        Uses scontrol and sinfo to get partition hardware specs.
        """
        try:
            # Get partition info via scontrol
            result = subprocess.run(
                ['scontrol', 'show', 'partition', partition],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0:
                logger.debug(f"scontrol failed for partition {partition}")
                return None
            
            partition_info = result.stdout
            
            # Parse TotalCPUs, TotalNodes to estimate CPUs per node
            # Note: This gives cluster totals, not per-node
            
            # Try sinfo for more specific node info
            result = subprocess.run(
                ['sinfo', '-p', partition, '-N', '-o', '%n %c %G %m', '--noheader'],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0 or not result.stdout.strip():
                # Try alternative format
                result = subprocess.run(
                    ['sinfo', '-p', partition, '-o', '%c %G %m', '--noheader'],
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
            gpu_vendor = GPUVendor.NONE
            memory_gb = 0.0
            
            for part in parts:
                # CPU cores (numeric value)
                if part.isdigit():
                    if cpu_cores is None:
                        cpu_cores = int(part)
                
                # GPU info (format: gpu:type:count or gpu:count)
                elif part.startswith('gpu:') or 'gpu' in part.lower():
                    gpu_match = re.search(r'gpu:(\w+:)?(\d+)', part, re.IGNORECASE)
                    if gpu_match:
                        gpus = int(gpu_match.group(2))
                        gpu_type = gpu_match.group(1) or ''
                        if 'a100' in gpu_type.lower() or 'v100' in gpu_type.lower():
                            gpu_vendor = GPUVendor.NVIDIA
                        elif 'mi' in gpu_type.lower():
                            gpu_vendor = GPUVendor.AMD
                
                # Memory (format: number or number+suffix)
                elif re.match(r'^\d+[MGT]?$', part):
                    mem_match = re.match(r'^(\d+)([MGT])?$', part)
                    if mem_match:
                        mem_val = int(mem_match.group(1))
                        suffix = mem_match.group(2) or 'M'
                        if suffix == 'G':
                            memory_gb = mem_val
                        elif suffix == 'T':
                            memory_gb = mem_val * 1024
                        else:
                            memory_gb = mem_val / 1024
            
            if cpu_cores is None:
                return None
            
            # If we couldn't detect GPU vendor, try to infer from partition name
            if gpus > 0 and gpu_vendor == GPUVendor.NONE:
                partition_lower = partition.lower()
                if 'nvidia' in partition_lower or 'cuda' in partition_lower:
                    gpu_vendor = GPUVendor.NVIDIA
                elif 'amd' in partition_lower or 'rocm' in partition_lower:
                    gpu_vendor = GPUVendor.AMD
                else:
                    # Default assumption for GPU partitions
                    gpu_vendor = GPUVendor.NVIDIA
            
            return NodeTopology(
                cpu_cores=cpu_cores,
                gpus=gpus,
                gpu_vendor=gpu_vendor,
                memory_gb=memory_gb,
            )
            
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
            logger.debug(f"Partition detection failed: {e}")
            return None
    
    def _detect_from_system(self) -> Optional[NodeTopology]:
        """
        Detect topology via direct system introspection.
        
        This is the fallback method that works without SLURM.
        """
        cpu_cores = self._detect_cpu_count()
        if cpu_cores is None:
            return None
        
        cpu_info = self._detect_cpu_topology()
        
        # GPU detection
        gpus = 0
        gpu_vendor = GPUVendor.NONE
        gpu_model = ""
        
        nvidia_info = self._detect_nvidia_gpus()
        if nvidia_info:
            gpus = nvidia_info['count']
            gpu_vendor = GPUVendor.NVIDIA
            gpu_model = nvidia_info.get('model', '')
        else:
            amd_info = self._detect_amd_gpus()
            if amd_info:
                gpus = amd_info['count']
                gpu_vendor = GPUVendor.AMD
                gpu_model = amd_info.get('model', '')
        
        memory_gb = self._detect_memory()
        
        return NodeTopology(
            cpu_cores=cpu_cores,
            cpu_sockets=cpu_info.get('sockets', 1),
            cores_per_socket=cpu_info.get('cores_per_socket', cpu_cores),
            threads_per_core=cpu_info.get('threads_per_core', 1),
            gpus=gpus,
            gpu_vendor=gpu_vendor,
            gpu_model=gpu_model,
            memory_gb=memory_gb,
        )
    
    def _detect_cpu_count(self) -> Optional[int]:
        """Detect CPU core count from system."""
        # Try /proc/cpuinfo
        try:
            with open('/proc/cpuinfo', 'r') as f:
                content = f.read()
                # Count processor entries
                count = len(re.findall(r'^processor\s*:', content, re.MULTILINE))
                if count > 0:
                    return count
        except (FileNotFoundError, PermissionError):
            pass
        
        # Try os.cpu_count()
        count = os.cpu_count()
        if count:
            return count
        
        # Try nproc
        try:
            result = subprocess.run(
                ['nproc', '--all'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                return int(result.stdout.strip())
        except (subprocess.TimeoutExpired, FileNotFoundError, ValueError):
            pass
        
        return None
    
    def _detect_cpu_topology(self) -> Dict[str, int]:
        """Detect detailed CPU topology."""
        info = {'sockets': 1, 'cores_per_socket': 1, 'threads_per_core': 1}
        
        # Try lscpu
        try:
            result = subprocess.run(
                ['lscpu'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                output = result.stdout
                
                # Parse lscpu output
                for line in output.split('\n'):
                    if ':' in line:
                        key, value = line.split(':', 1)
                        key = key.strip().lower()
                        value = value.strip()
                        
                        if 'socket' in key and value.isdigit():
                            info['sockets'] = int(value)
                        elif 'core(s) per socket' in key and value.isdigit():
                            info['cores_per_socket'] = int(value)
                        elif 'thread(s) per core' in key and value.isdigit():
                            info['threads_per_core'] = int(value)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        return info
    
    def _detect_nvidia_gpus(self) -> Optional[Dict]:
        """Detect NVIDIA GPUs using nvidia-smi."""
        try:
            result = subprocess.run(
                ['nvidia-smi', '--query-gpu=count,name,memory.total', 
                 '--format=csv,noheader,nounits'],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0 and result.stdout.strip():
                lines = result.stdout.strip().split('\n')
                # First line has count for first GPU, count lines for total
                count = len(lines)
                
                # Parse first line for model info
                parts = lines[0].split(',')
                model = parts[1].strip() if len(parts) > 1 else ''
                memory = 0
                if len(parts) > 2:
                    try:
                        memory = float(parts[2].strip()) / 1024  # Convert MB to GB
                    except ValueError:
                        pass
                
                return {
                    'count': count,
                    'model': model,
                    'memory_gb': memory
                }
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        # Alternative: check /dev/nvidia*
        try:
            import glob
            nvidia_devs = glob.glob('/dev/nvidia[0-9]*')
            if nvidia_devs:
                return {'count': len(nvidia_devs), 'model': '', 'memory_gb': 0}
        except Exception:
            pass
        
        return None
    
    def _detect_amd_gpus(self) -> Optional[Dict]:
        """Detect AMD GPUs using rocm-smi."""
        try:
            result = subprocess.run(
                ['rocm-smi', '--showproductname'],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0:
                # Count GPU entries
                count = len(re.findall(r'GPU\[(\d+)\]', result.stdout))
                if count > 0:
                    # Try to get model
                    model_match = re.search(r'Card series:\s*(.+)', result.stdout)
                    model = model_match.group(1).strip() if model_match else ''
                    return {'count': count, 'model': model, 'memory_gb': 0}
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        # Alternative: check /dev/dri/renderD*
        try:
            import glob
            render_devs = glob.glob('/dev/dri/renderD*')
            if render_devs:
                # This is a rough heuristic - not all render devices are GPUs
                return {'count': len(render_devs), 'model': '', 'memory_gb': 0}
        except Exception:
            pass
        
        return None
    
    def _detect_memory(self) -> float:
        """Detect system memory in GB."""
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemTotal:'):
                        # Format: "MemTotal:       123456789 kB"
                        parts = line.split()
                        if len(parts) >= 2:
                            kb = int(parts[1])
                            return kb / (1024 * 1024)  # Convert to GB
        except (FileNotFoundError, PermissionError, ValueError):
            pass
        
        return 0.0
    
    def compute_mpi_mapping(
        self,
        topology: NodeTopology,
        num_nodes: int = 1,
        user_ranks_per_node: Optional[int] = None,
        user_cores_per_rank: Optional[int] = None,
    ) -> MPIMapping:
        """
        Compute optimal MPI mapping based on topology.
        
        This implements the core logic for deriving MPI configuration from
        hardware topology. The algorithm prioritizes GPU utilization when
        GPUs are present.
        
        Algorithm:
        1. If GPUs present: ranks_per_node = gpus (1 rank per GPU)
        2. cores_per_rank = cpu_cores / ranks_per_node
        3. Validate that mapping fits within hardware constraints
        
        Args:
            topology: Detected hardware topology
            num_nodes: Number of nodes in the job
            user_ranks_per_node: User override for ranks per node
            user_cores_per_rank: User override for cores per rank
        
        Returns:
            MPIMapping with computed configuration
        
        Raises:
            ValueError: If configuration is invalid for the topology
        """
        # Determine ranks per node
        if user_ranks_per_node is not None:
            ranks_per_node = user_ranks_per_node
        elif topology.gpus > 0:
            # GPU job: 1 rank per GPU (optimal for most GPU codes)
            ranks_per_node = topology.gpus
        else:
            # CPU job: 1 rank per physical core
            ranks_per_node = topology.physical_cores or topology.cpu_cores
        
        # Determine cores per rank
        if user_cores_per_rank is not None:
            cores_per_rank = user_cores_per_rank
        else:
            # Distribute cores evenly among ranks
            cores_per_rank = topology.cpu_cores // ranks_per_node
        
        # Determine GPUs per rank
        if topology.gpus > 0:
            gpus_per_rank = topology.gpus // ranks_per_node
        else:
            gpus_per_rank = 0
        
        # Total ranks
        total_ranks = ranks_per_node * num_nodes
        
        mapping = MPIMapping(
            ranks_per_node=ranks_per_node,
            cores_per_rank=cores_per_rank,
            gpus_per_rank=gpus_per_rank,
            total_ranks=total_ranks,
            binding_strategy="core" if cores_per_rank > 1 else "none"
        )
        
        # Validate
        is_valid, error = mapping.validate(topology)
        if not is_valid:
            raise ValueError(f"Invalid MPI mapping: {error}")
        
        logger.info(f"Computed MPI mapping for {num_nodes} node(s):")
        logger.info(f"  Ranks per node: {mapping.ranks_per_node}")
        logger.info(f"  Cores per rank: {mapping.cores_per_rank}")
        logger.info(f"  GPUs per rank: {mapping.gpus_per_rank}")
        logger.info(f"  Total ranks: {mapping.total_ranks}")
        
        return mapping


# Global detector instance (lazy initialization)
_global_detector: Optional[TopologyDetector] = None


def get_topology_detector() -> TopologyDetector:
    """Get the global topology detector instance."""
    global _global_detector
    if _global_detector is None:
        _global_detector = TopologyDetector()
    return _global_detector


def detect_topology(partition: Optional[str] = None) -> NodeTopology:
    """
    Convenience function to detect topology.
    
    Args:
        partition: Optional SLURM partition name
    
    Returns:
        NodeTopology with detected hardware configuration
    """
    return get_topology_detector().detect(partition)


def compute_mpi_mapping(
    num_nodes: int = 1,
    partition: Optional[str] = None,
    ranks_per_node: Optional[int] = None,
    cores_per_rank: Optional[int] = None,
) -> Tuple[NodeTopology, MPIMapping]:
    """
    Convenience function to detect topology and compute MPI mapping.
    
    Args:
        num_nodes: Number of nodes
        partition: Optional SLURM partition name
        ranks_per_node: Optional user override
        cores_per_rank: Optional user override
    
    Returns:
        Tuple of (NodeTopology, MPIMapping)
    """
    detector = get_topology_detector()
    topology = detector.detect(partition)
    mapping = detector.compute_mpi_mapping(
        topology,
        num_nodes=num_nodes,
        user_ranks_per_node=ranks_per_node,
        user_cores_per_rank=cores_per_rank,
    )
    return topology, mapping
