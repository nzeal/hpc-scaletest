#!/usr/bin/env python3
"""
GPU Detection and Configuration Module

Automatically detects GPU hardware and configures MPI tasks accordingly.
Supports NVIDIA, AMD, Intel GPUs across different HPC systems.
"""

import subprocess
import logging
import os
from typing import Dict, List, Optional, Tuple
from pathlib import Path
from dataclasses import dataclass


logger = logging.getLogger(__name__)


@dataclass
class GPUInfo:
    """Information about GPU hardware."""
    vendor: str  # 'nvidia', 'amd', 'intel', 'unknown'
    count_per_node: int
    model: str
    architecture: str  # 'cuda', 'rocm', 'oneapi', 'unknown'
    device_binding_var: str  # Environment variable for device binding
    total_memory_gb: Optional[float] = None


class GPUDetector:
    """Detect GPU hardware and provide configuration."""
    
    def __init__(self):
        self.gpu_info: Optional[GPUInfo] = None
    
    def detect(self, partition_name: Optional[str] = None) -> Optional[GPUInfo]:
        """
        Detect GPU hardware on the system.
        
        Args:
            partition_name: Optional partition name to query partition-specific info
            
        Returns:
            GPUInfo object if GPUs detected, None otherwise
        """
        # Try detection methods in order of reliability
        methods = [
            self._detect_from_slurm_partition,
            self._detect_nvidia_smi,
            self._detect_rocm_smi,
            self._detect_intel_gpu,
            self._detect_from_environment,
        ]
        
        for method in methods:
            try:
                if partition_name and method == self._detect_from_slurm_partition:
                    gpu_info = method(partition_name)
                else:
                    gpu_info = method()
                    
                if gpu_info:
                    logger.info(f"✓ Detected {gpu_info.count_per_node} {gpu_info.vendor.upper()} "
                              f"GPU(s) per node ({gpu_info.model})")
                    self.gpu_info = gpu_info
                    return gpu_info
            except Exception as e:
                logger.debug(f"Detection method {method.__name__} failed: {e}")
                continue
        
        logger.debug("No GPUs detected")
        return None
    
    def _detect_from_slurm_partition(self, partition_name: str) -> Optional[GPUInfo]:
        """Detect GPU info from SLURM partition configuration."""
        try:
            # Query partition info with scontrol
            result = subprocess.run(
                ['scontrol', 'show', 'partition', partition_name],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode != 0:
                return None
            
            # Parse GRES (Generic RESources) field
            for line in result.stdout.split('\n'):
                if 'TRES=' in line or 'Gres=' in line:
                    # Look for gpu specifications like "gpu:4" or "gpu:a100:4"
                    import re
                    gpu_match = re.search(r'gpu(?::(\w+))?:(\d+)', line, re.IGNORECASE)
                    if gpu_match:
                        model = gpu_match.group(1) or 'unknown'
                        count = int(gpu_match.group(2))
                        
                        # Determine vendor from model name
                        vendor = 'nvidia'  # Default
                        if any(x in model.lower() for x in ['mi', 'gfx', 'vega', 'radeon']):
                            vendor = 'amd'
                        elif any(x in model.lower() for x in ['pvc', 'ponte', 'xe']):
                            vendor = 'intel'
                        
                        return GPUInfo(
                            vendor=vendor,
                            count_per_node=count,
                            model=model,
                            architecture=self._get_architecture(vendor),
                            device_binding_var=self._get_binding_var(vendor)
                        )
        except Exception as e:
            logger.debug(f"SLURM partition detection failed: {e}")
        
        return None
    
    def _detect_nvidia_smi(self) -> Optional[GPUInfo]:
        """Detect NVIDIA GPUs using nvidia-smi."""
        try:
            result = subprocess.run(
                ['nvidia-smi', '--query-gpu=count,name,memory.total', '--format=csv,noheader'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0 and result.stdout.strip():
                lines = result.stdout.strip().split('\n')
                if lines:
                    # Parse first GPU info
                    parts = lines[0].split(',')
                    count = len(lines)  # Number of GPUs
                    model = parts[1].strip() if len(parts) > 1 else 'NVIDIA GPU'
                    
                    # Get memory if available
                    memory_gb = None
                    if len(parts) > 2:
                        try:
                            memory_str = parts[2].strip().split()[0]
                            memory_gb = float(memory_str) / 1024  # Convert MB to GB
                        except:
                            pass
                    
                    return GPUInfo(
                        vendor='nvidia',
                        count_per_node=count,
                        model=model,
                        architecture='cuda',
                        device_binding_var='CUDA_VISIBLE_DEVICES',
                        total_memory_gb=memory_gb
                    )
        except FileNotFoundError:
            logger.debug("nvidia-smi not found")
        except Exception as e:
            logger.debug(f"nvidia-smi detection failed: {e}")
        
        return None
    
    def _detect_rocm_smi(self) -> Optional[GPUInfo]:
        """Detect AMD GPUs using rocm-smi."""
        try:
            result = subprocess.run(
                ['rocm-smi', '--showproductname'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0 and result.stdout.strip():
                # Count GPU entries
                lines = [l for l in result.stdout.split('\n') if 'GPU' in l]
                count = len(lines)
                
                # Try to get model name
                model = 'AMD GPU'
                for line in lines:
                    if 'Card series:' in line or 'Card model:' in line:
                        model = line.split(':')[-1].strip()
                        break
                
                return GPUInfo(
                    vendor='amd',
                    count_per_node=count,
                    model=model,
                    architecture='rocm',
                    device_binding_var='ROCR_VISIBLE_DEVICES'
                )
        except FileNotFoundError:
            logger.debug("rocm-smi not found")
        except Exception as e:
            logger.debug(f"rocm-smi detection failed: {e}")
        
        return None
    
    def _detect_intel_gpu(self) -> Optional[GPUInfo]:
        """Detect Intel GPUs."""
        try:
            # Try intel_gpu_top or similar tools
            result = subprocess.run(
                ['clinfo'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0:
                # Look for Intel GPU devices
                intel_count = result.stdout.count('Intel') if 'Intel' in result.stdout else 0
                if intel_count > 0:
                    return GPUInfo(
                        vendor='intel',
                        count_per_node=intel_count,
                        model='Intel GPU',
                        architecture='oneapi',
                        device_binding_var='ZE_AFFINITY_MASK'
                    )
        except FileNotFoundError:
            pass
        except Exception as e:
            logger.debug(f"Intel GPU detection failed: {e}")
        
        return None
    
    def _detect_from_environment(self) -> Optional[GPUInfo]:
        """Detect GPUs from environment variables."""
        # Check common GPU-related environment variables
        gpu_vars = {
            'CUDA_VISIBLE_DEVICES': ('nvidia', 'cuda'),
            'ROCR_VISIBLE_DEVICES': ('amd', 'rocm'),
            'GPU_DEVICE_ORDINAL': ('nvidia', 'cuda'),
        }
        
        for var, (vendor, arch) in gpu_vars.items():
            if var in os.environ:
                devices = os.environ[var]
                if devices and devices != '-1':
                    count = len(devices.split(','))
                    return GPUInfo(
                        vendor=vendor,
                        count_per_node=count,
                        model=f'{vendor.upper()} GPU',
                        architecture=arch,
                        device_binding_var=var
                    )
        
        return None
    
    def _get_architecture(self, vendor: str) -> str:
        """Get GPU architecture from vendor."""
        arch_map = {
            'nvidia': 'cuda',
            'amd': 'rocm',
            'intel': 'oneapi'
        }
        return arch_map.get(vendor.lower(), 'unknown')
    
    def _get_binding_var(self, vendor: str) -> str:
        """Get device binding environment variable for vendor."""
        binding_vars = {
            'nvidia': 'CUDA_VISIBLE_DEVICES',
            'amd': 'ROCR_VISIBLE_DEVICES',
            'intel': 'ZE_AFFINITY_MASK'
        }
        return binding_vars.get(vendor.lower(), 'GPU_DEVICE_ORDINAL')
    
    def get_optimal_tasks_per_node(self, 
                                   gpu_info: Optional[GPUInfo] = None,
                                   requested_tasks: Optional[int] = None) -> int:
        """
        Determine optimal number of MPI tasks per node.
        
        For GPU systems: matches GPU count (1 task per GPU)
        For CPU systems: uses requested value or system default
        
        Args:
            gpu_info: GPU information (uses cached if None)
            requested_tasks: User-requested tasks per node
            
        Returns:
            Optimal number of tasks per node
        """
        if gpu_info is None:
            gpu_info = self.gpu_info
        
        if gpu_info is not None:
            # For GPU jobs: 1 MPI task per GPU
            optimal = gpu_info.count_per_node
            
            if requested_tasks and requested_tasks != optimal:
                logger.warning(
                    f"⚠ Requested {requested_tasks} tasks/node but system has "
                    f"{optimal} GPUs. Using {optimal} for optimal GPU binding."
                )
            
            return optimal
        else:
            # For CPU jobs: use requested or default
            return requested_tasks if requested_tasks else 128
    
    def generate_gpu_binding(self, 
                            rank: int,
                            tasks_per_node: int,
                            gpu_info: Optional[GPUInfo] = None) -> Dict[str, str]:
        """
        Generate GPU binding environment variables for a given MPI rank.
        
        Args:
            rank: MPI rank number
            tasks_per_node: Number of tasks per node
            gpu_info: GPU information
            
        Returns:
            Dictionary of environment variables for GPU binding
        """
        if gpu_info is None:
            gpu_info = self.gpu_info
        
        if gpu_info is None:
            return {}
        
        # Calculate which GPU this rank should use
        local_rank = rank % tasks_per_node
        gpu_id = local_rank % gpu_info.count_per_node
        
        env_vars = {
            gpu_info.device_binding_var: str(gpu_id)
        }
        
        # Add vendor-specific variables
        if gpu_info.vendor == 'nvidia':
            env_vars['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'
        elif gpu_info.vendor == 'amd':
            env_vars['HIP_VISIBLE_DEVICES'] = str(gpu_id)
        
        return env_vars
    
    def validate_gpu_configuration(self,
                                  tasks_per_node: int,
                                  nodes: int,
                                  gpu_info: Optional[GPUInfo] = None) -> Tuple[bool, str]:
        """
        Validate that the configuration is suitable for GPU execution.
        
        Returns:
            (is_valid, message)
        """
        if gpu_info is None:
            gpu_info = self.gpu_info
        
        if gpu_info is None:
            return True, "No GPU configuration to validate"
        
        # Check if tasks_per_node matches GPU count
        if tasks_per_node != gpu_info.count_per_node:
            return False, (
                f"Mismatch: {tasks_per_node} tasks/node but {gpu_info.count_per_node} GPUs/node. "
                f"For optimal GPU utilization, use {gpu_info.count_per_node} tasks/node."
            )
        
        # Check if total GPU count is reasonable
        total_gpus = tasks_per_node * nodes
        if total_gpus > 1024:  # Arbitrary large number check
            return False, f"Configuration would use {total_gpus} GPUs which seems excessive"
        
        return True, f"✓ Configuration validated: {total_gpus} GPUs across {nodes} nodes"


def get_gpu_info_for_partition(partition_name: str) -> Optional[GPUInfo]:
    """
    Convenience function to detect GPU info for a specific partition.
    
    Args:
        partition_name: Name of the SLURM partition
        
    Returns:
        GPUInfo if GPUs detected, None otherwise
    """
    detector = GPUDetector()
    return detector.detect(partition_name=partition_name)


def get_optimal_gpu_tasks(partition_name: Optional[str] = None,
                         requested_tasks: Optional[int] = None) -> int:
    """
    Get optimal number of tasks per node for GPU execution.
    
    Args:
        partition_name: Optional partition name
        requested_tasks: User-requested tasks
        
    Returns:
        Optimal number of tasks per node
    """
    detector = GPUDetector()
    gpu_info = detector.detect(partition_name=partition_name)
    return detector.get_optimal_tasks_per_node(gpu_info, requested_tasks)


if __name__ == '__main__':
    # Test GPU detection
    logging.basicConfig(level=logging.INFO)
    
    detector = GPUDetector()
    gpu_info = detector.detect()
    
    if gpu_info:
        print(f"\n{'='*60}")
        print("GPU Detection Results:")
        print(f"{'='*60}")
        print(f"Vendor:           {gpu_info.vendor.upper()}")
        print(f"Model:            {gpu_info.model}")
        print(f"GPUs per node:    {gpu_info.count_per_node}")
        print(f"Architecture:     {gpu_info.architecture}")
        print(f"Binding variable: {gpu_info.device_binding_var}")
        if gpu_info.total_memory_gb:
            print(f"Memory per GPU:   {gpu_info.total_memory_gb:.1f} GB")
        print(f"{'='*60}")
        
        # Show optimal configuration
        optimal_tasks = detector.get_optimal_tasks_per_node()
        print(f"\nOptimal MPI tasks per node: {optimal_tasks}")
        
        # Show example binding
        print(f"\nExample GPU binding for first 3 ranks:")
        for rank in range(3):
            binding = detector.generate_gpu_binding(rank, optimal_tasks)
            print(f"  Rank {rank}: {binding}")
    else:
        print("No GPUs detected on this system")
