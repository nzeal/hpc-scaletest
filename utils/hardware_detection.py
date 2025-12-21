#!/usr/bin/env python3
"""
Hardware-Aware Configuration Manager

Automatically configures HPC jobs based on detected hardware:
- GPU count and type (NVIDIA, AMD, Intel)
- CPU architecture (x86_64, ARM, RISC-V)
- Node topology and memory

Ensures optimal resource utilization across different systems.
"""

import logging
import subprocess
import platform
from typing import Dict, Optional, Tuple, Any
from dataclasses import dataclass, field

from utils.gpu_detection import GPUDetector, GPUInfo
from utils.partition_validator import PartitionValidator


logger = logging.getLogger(__name__)


@dataclass
class HardwareConfig:
    """Complete hardware configuration for a node."""
    # Architecture
    architecture: str  # 'x86_64', 'aarch64', 'riscv64'
    cpu_vendor: str    # 'intel', 'amd', 'arm', 'riscv', 'unknown'
    
    # CPU resources
    cpus_per_node: int
    cores_per_socket: int
    sockets_per_node: int
    threads_per_core: int
    
    # GPU resources
    has_gpu: bool = False
    gpu_info: Optional[GPUInfo] = None
    
    # Memory
    memory_per_node_gb: Optional[float] = None
    
    # Network
    network_fabric: Optional[str] = None  # 'infiniband', 'ethernet', 'omnipath'
    
    # Derived properties
    optimal_tasks_per_node: int = field(init=False)
    
    def __post_init__(self):
        """Calculate derived properties."""
        if self.has_gpu and self.gpu_info:
            # For GPU: 1 task per GPU
            self.optimal_tasks_per_node = self.gpu_info.count_per_node
        else:
            # For CPU: use physical cores (not hyperthreads)
            self.optimal_tasks_per_node = self.cores_per_socket * self.sockets_per_node


class HardwareDetector:
    """Detect hardware configuration."""
    
    def __init__(self):
        self.gpu_detector = GPUDetector()
        self.partition_validator = PartitionValidator()
    
    def detect_full_config(self, 
                          partition_name: Optional[str] = None,
                          hardware_type: str = 'cpu') -> HardwareConfig:
        """
        Detect complete hardware configuration.
        
        Args:
            partition_name: Optional partition name for partition-specific detection
            hardware_type: 'cpu' or 'gpu'
            
        Returns:
            HardwareConfig object with detected hardware
        """
        # Detect CPU architecture
        arch = platform.machine()
        cpu_vendor = self._detect_cpu_vendor()
        
        # Detect CPU topology
        cpus, cores, sockets, threads = self._detect_cpu_topology()
        
        # Detect GPU if requested
        gpu_info = None
        has_gpu = False
        if hardware_type == 'gpu':
            gpu_info = self.gpu_detector.detect(partition_name=partition_name)
            has_gpu = gpu_info is not None
        
        # Detect memory
        memory_gb = self._detect_memory()
        
        # Detect network
        network = self._detect_network()
        
        config = HardwareConfig(
            architecture=arch,
            cpu_vendor=cpu_vendor,
            cpus_per_node=cpus,
            cores_per_socket=cores,
            sockets_per_node=sockets,
            threads_per_core=threads,
            has_gpu=has_gpu,
            gpu_info=gpu_info,
            memory_per_node_gb=memory_gb,
            network_fabric=network
        )
        
        logger.info(f"✓ Hardware detected: {config.architecture} "
                   f"({config.cpu_vendor}, {config.cpus_per_node} CPUs"
                   f"{f', {config.gpu_info.count_per_node} GPUs' if has_gpu else ''})")
        
        return config
    
    def _detect_cpu_vendor(self) -> str:
        """Detect CPU vendor."""
        try:
            # Try lscpu first
            result = subprocess.run(
                ['lscpu'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0:
                output = result.stdout.lower()
                if 'genuineintel' in output or 'intel' in output:
                    return 'intel'
                elif 'authenticamd' in output or 'amd' in output:
                    return 'amd'
                elif 'arm' in output:
                    return 'arm'
                elif 'riscv' in output:
                    return 'riscv'
        except:
            pass
        
        # Fallback to platform
        machine = platform.machine().lower()
        if 'arm' in machine or 'aarch' in machine:
            return 'arm'
        elif 'riscv' in machine:
            return 'riscv'
        
        return 'unknown'
    
    def _detect_cpu_topology(self) -> Tuple[int, int, int, int]:
        """
        Detect CPU topology.
        
        Returns:
            (cpus_per_node, cores_per_socket, sockets_per_node, threads_per_core)
        """
        cpus = cores = sockets = threads = 1
        
        try:
            # Try lscpu
            result = subprocess.run(
                ['lscpu'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0:
                for line in result.stdout.split('\n'):
                    if 'CPU(s):' in line and 'NUMA' not in line and 'On-line' not in line:
                        cpus = int(line.split(':')[1].strip())
                    elif 'Core(s) per socket:' in line:
                        cores = int(line.split(':')[1].strip())
                    elif 'Socket(s):' in line:
                        sockets = int(line.split(':')[1].strip())
                    elif 'Thread(s) per core:' in line:
                        threads = int(line.split(':')[1].strip())
        except:
            # Fallback to simple detection
            try:
                import os
                cpus = os.cpu_count() or 1
            except:
                cpus = 1
        
        # Ensure logical consistency
        if cores == 1 and sockets == 1:
            # Estimate from total CPUs
            cores = cpus // threads if threads > 1 else cpus
            sockets = 1
        
        return cpus, cores, sockets, threads
    
    def _detect_memory(self) -> Optional[float]:
        """Detect total memory in GB."""
        try:
            result = subprocess.run(
                ['free', '-g'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0:
                lines = result.stdout.split('\n')
                for line in lines:
                    if 'Mem:' in line:
                        parts = line.split()
                        if len(parts) >= 2:
                            return float(parts[1])
        except:
            pass
        
        return None
    
    def _detect_network(self) -> Optional[str]:
        """Detect network fabric."""
        try:
            # Check for InfiniBand
            result = subprocess.run(
                ['ibstat'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0 and 'Infiniband' in result.stdout:
                return 'infiniband'
        except:
            pass
        
        try:
            # Check for OmniPath
            result = subprocess.run(
                ['opainfo'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                return 'omnipath'
        except:
            pass
        
        return 'ethernet'  # Default assumption
    
    def get_optimal_mpi_config(self,
                              config: HardwareConfig,
                              nodes: int,
                              scaling_type: str = 'strong') -> Dict[str, Any]:
        """
        Get optimal MPI configuration for the hardware.
        
        Args:
            config: Hardware configuration
            nodes: Number of nodes
            scaling_type: 'strong' or 'weak'
            
        Returns:
            Dictionary with MPI configuration
        """
        tasks_per_node = config.optimal_tasks_per_node
        total_tasks = tasks_per_node * nodes
        
        # Calculate process topology
        if config.has_gpu and config.gpu_info:
            # GPU configuration: match GPU layout
            topology = self._calculate_gpu_topology(config.gpu_info.count_per_node, nodes)
        else:
            # CPU configuration: optimize for cache/NUMA
            topology = self._calculate_cpu_topology(tasks_per_node, nodes)
        
        mpi_config = {
            'tasks_per_node': tasks_per_node,
            'total_tasks': total_tasks,
            'nodes': nodes,
            'topology': topology,
            'binding_strategy': self._get_binding_strategy(config),
            'environment': self._get_mpi_environment(config)
        }
        
        return mpi_config
    
    def _calculate_gpu_topology(self, gpus_per_node: int, nodes: int) -> Tuple[int, int, int]:
        """
        Calculate process topology for GPU layout.
        
        For GPUs, we want:
        - Each dimension to be a factor that distributes GPUs evenly
        - Prefer 2D layouts for small GPU counts
        - Use 3D for larger configurations
        """
        total_gpus = gpus_per_node * nodes
        
        # Simple factorization
        if total_gpus <= 4:
            # 1D or 2D layout
            return (total_gpus, 1, 1)
        elif total_gpus <= 16:
            # 2D layout
            import math
            sqrt = int(math.sqrt(total_gpus))
            return (sqrt, total_gpus // sqrt, 1)
        else:
            # 3D layout
            return self._factorize_3d(total_gpus)
    
    def _calculate_cpu_topology(self, tasks_per_node: int, nodes: int) -> Tuple[int, int, int]:
        """Calculate process topology for CPU layout."""
        total_tasks = tasks_per_node * nodes
        return self._factorize_3d(total_tasks)
    
    def _factorize_3d(self, n: int) -> Tuple[int, int, int]:
        """Factorize n into three roughly equal factors."""
        import math
        
        # Try to find three factors close to cube root
        target = int(math.pow(n, 1/3))
        
        for z in range(target, 0, -1):
            if n % z == 0:
                remaining = n // z
                sqrt = int(math.sqrt(remaining))
                for y in range(sqrt, 0, -1):
                    if remaining % y == 0:
                        x = remaining // y
                        return (x, y, z)
        
        # Fallback: unbalanced factorization
        return (n, 1, 1)
    
    def _get_binding_strategy(self, config: HardwareConfig) -> str:
        """Get optimal process binding strategy."""
        if config.has_gpu:
            return 'gpu'  # Bind to GPUs
        elif config.threads_per_core > 1:
            return 'core'  # Bind to physical cores
        else:
            return 'socket'  # Bind to sockets
    
    def _get_mpi_environment(self, config: HardwareConfig) -> Dict[str, str]:
        """Get MPI environment variables for optimal performance."""
        env = {}
        
        # Common MPI settings
        if config.threads_per_core > 1:
            env['OMP_NUM_THREADS'] = '1'  # Disable OpenMP if using all cores
        
        # GPU-specific settings
        if config.has_gpu and config.gpu_info:
            if config.gpu_info.vendor == 'nvidia':
                env['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'
            elif config.gpu_info.vendor == 'amd':
                env['HIP_LAUNCH_BLOCKING'] = '0'
        
        # Network-specific settings
        if config.network_fabric == 'infiniband':
            env['UCX_TLS'] = 'rc,ud,sm,self'
            env['UCX_NET_DEVICES'] = 'mlx5_0:1'
        
        return env


def get_hardware_config(partition_name: Optional[str] = None,
                       hardware_type: str = 'cpu') -> HardwareConfig:
    """
    Convenience function to detect hardware configuration.
    
    Args:
        partition_name: Optional partition name
        hardware_type: 'cpu' or 'gpu'
        
    Returns:
        HardwareConfig object
    """
    detector = HardwareDetector()
    return detector.detect_full_config(partition_name, hardware_type)


if __name__ == '__main__':
    # Test hardware detection
    logging.basicConfig(level=logging.INFO)
    
    detector = HardwareDetector()
    
    # Detect CPU configuration
    print("\nCPU Configuration:")
    print("="*60)
    cpu_config = detector.detect_full_config(hardware_type='cpu')
    print(f"Architecture:     {cpu_config.architecture}")
    print(f"CPU Vendor:       {cpu_config.cpu_vendor}")
    print(f"CPUs per node:    {cpu_config.cpus_per_node}")
    print(f"Cores per socket: {cpu_config.cores_per_socket}")
    print(f"Sockets:          {cpu_config.sockets_per_node}")
    print(f"Threads per core: {cpu_config.threads_per_core}")
    print(f"Optimal tasks:    {cpu_config.optimal_tasks_per_node}")
    if cpu_config.memory_per_node_gb:
        print(f"Memory:           {cpu_config.memory_per_node_gb:.0f} GB")
    if cpu_config.network_fabric:
        print(f"Network:          {cpu_config.network_fabric}")
    
    # Detect GPU configuration
    print("\nGPU Configuration:")
    print("="*60)
    gpu_config = detector.detect_full_config(hardware_type='gpu')
    if gpu_config.has_gpu:
        print(f"GPUs detected:    Yes")
        print(f"GPU vendor:       {gpu_config.gpu_info.vendor.upper()}")
        print(f"GPU model:        {gpu_config.gpu_info.model}")
        print(f"GPUs per node:    {gpu_config.gpu_info.count_per_node}")
        print(f"Architecture:     {gpu_config.gpu_info.architecture}")
        print(f"Optimal tasks:    {gpu_config.optimal_tasks_per_node}")
    else:
        print("No GPUs detected")
    
    print("="*60)
