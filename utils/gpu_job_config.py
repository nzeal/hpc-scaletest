#!/usr/bin/env python3
"""
GPU-Aware Job Configuration

Integrates GPU detection into job configuration workflow.
Ensures correct SLURM directives for GPU jobs on any system.
"""

import logging
from typing import Dict, Any, Optional, Tuple
from dataclasses import dataclass

from utils.gpu_detection import GPUDetector, GPUInfo
from utils.partition_validator import PartitionValidator


logger = logging.getLogger(__name__)


@dataclass
class GPUJobConfig:
    """GPU-specific job configuration."""
    gpus_per_node: int
    tasks_per_node: int
    gpu_binding: str = "closest"
    gpu_vendor: str = "nvidia"
    device_var: str = "CUDA_VISIBLE_DEVICES"


class GPUConfigManager:
    """Manages GPU configuration for jobs."""
    
    def __init__(self):
        self.gpu_detector = GPUDetector()
        self.partition_validator = PartitionValidator()
        self.gpu_info: Optional[GPUInfo] = None
    
    def configure_for_partition(
        self,
        partition: str,
        hardware_type: str = 'cpu',
        requested_tasks_per_node: Optional[int] = None
    ) -> Tuple[int, int]:
        """
        Configure GPU/CPU tasks for a given partition.
        
        Args:
            partition: SLURM partition name
            hardware_type: 'cpu' or 'gpu'
            requested_tasks_per_node: User-requested tasks (optional)
            
        Returns:
            (tasks_per_node, gpus_per_node)
        """
        # Validate partition first
        is_valid, message = self.partition_validator.validate_partition(partition)
        if not is_valid:
            logger.error(message)
            raise ValueError(f"Invalid partition: {partition}")
        
        logger.info(f"✓ Partition '{partition}' validated")
        
        # For CPU jobs, use requested or default
        if hardware_type == 'cpu':
            tasks_per_node = requested_tasks_per_node if requested_tasks_per_node else 128
            return tasks_per_node, 0
        
        # For GPU jobs, detect hardware
        self.gpu_info = self.gpu_detector.detect(partition_name=partition)
        
        if not self.gpu_info:
            logger.warning(f"⚠ Could not detect GPUs for partition '{partition}'")
            logger.warning("   Falling back to user-specified configuration")
            tasks_per_node = requested_tasks_per_node if requested_tasks_per_node else 4
            gpus_per_node = requested_tasks_per_node if requested_tasks_per_node else 4
            return tasks_per_node, gpus_per_node
        
        # GPU detected - use 1 task per GPU
        gpus_per_node = self.gpu_info.count_per_node
        tasks_per_node = gpus_per_node  # 1:1 mapping
        
        logger.info(f"✓ Detected {gpus_per_node} {self.gpu_info.vendor.upper()} "
                   f"GPU(s) per node ({self.gpu_info.model})")
        logger.info(f"✓ Setting tasks_per_node = {tasks_per_node} (1 task per GPU)")
        
        # Warn if user requested different
        if requested_tasks_per_node and requested_tasks_per_node != tasks_per_node:
            logger.warning(f"⚠ User requested {requested_tasks_per_node} tasks/node, "
                         f"but GPU count is {gpus_per_node}")
            logger.warning(f"   Using {tasks_per_node} tasks/node for optimal GPU utilization")
        
        return tasks_per_node, gpus_per_node
    
    def get_gpu_job_config(self) -> Optional[GPUJobConfig]:
        """Get GPU-specific configuration."""
        if not self.gpu_info:
            return None
        
        return GPUJobConfig(
            gpus_per_node=self.gpu_info.count_per_node,
            tasks_per_node=self.gpu_info.count_per_node,
            gpu_binding="closest",
            gpu_vendor=self.gpu_info.vendor,
            device_var=self.gpu_info.device_binding_var
        )
    
    def generate_gpu_env_vars(self, rank: int, tasks_per_node: int) -> Dict[str, str]:
        """Generate GPU environment variables for a given rank."""
        if not self.gpu_info:
            return {}
        
        return self.gpu_detector.generate_gpu_binding(
            rank, tasks_per_node, self.gpu_info
        )
    
    def validate_scaling_config(
        self,
        initial_procs: list,
        scaling_type: str,
        nodes: int
    ) -> Tuple[bool, str]:
        """
        Validate scaling configuration for GPUs.
        
        Args:
            initial_procs: Initial MPI decomposition [x, y, z]
            scaling_type: 'weak' or 'strong'
            nodes: Number of nodes
            
        Returns:
            (is_valid, message)
        """
        if not self.gpu_info:
            return True, "No GPU validation needed"
        
        # Calculate initial task count
        initial_tasks = 1
        for p in initial_procs:
            initial_tasks *= p
        
        tasks_per_node = self.gpu_info.count_per_node
        
        # For weak scaling, initial tasks should match GPU count
        if scaling_type == 'weak':
            if initial_tasks != tasks_per_node:
                return False, (
                    f"Weak scaling: initial task count ({initial_tasks}) "
                    f"must match GPU count ({tasks_per_node}) per node.\n"
                    f"Suggested initial_procs for {tasks_per_node} GPUs:\n"
                    f"  - 2D: [{tasks_per_node}, 1, 1] or [{tasks_per_node//2}, 2, 1]\n"
                    f"  - 3D: Factorize {tasks_per_node} into 3 dimensions"
                )
        
        # For strong scaling, verify divisibility
        if scaling_type == 'strong':
            total_gpus = tasks_per_node * nodes
            if initial_tasks > total_gpus:
                return False, (
                    f"Strong scaling: initial task count ({initial_tasks}) "
                    f"exceeds total GPU count ({total_gpus})"
                )
        
        return True, f"✓ Scaling configuration valid for {tasks_per_node} GPUs/node"


def configure_gpu_job(
    partition: str,
    hardware_type: str = 'gpu',
    procs_per_node: Optional[int] = None,
    initial_procs: Optional[list] = None,
    scaling_type: str = 'weak',
    nodes: int = 1
) -> Dict[str, Any]:
    """
    Configure GPU job with automatic hardware detection.
    
    Args:
        partition: SLURM partition name
        hardware_type: 'cpu' or 'gpu'
        procs_per_node: User-requested tasks per node (optional)
        initial_procs: Initial MPI decomposition for scaling
        scaling_type: 'weak' or 'strong'
        nodes: Number of nodes
        
    Returns:
        Configuration dictionary with GPU settings
    """
    manager = GPUConfigManager()
    
    # Configure based on partition
    tasks_per_node, gpus_per_node = manager.configure_for_partition(
        partition, hardware_type, procs_per_node
    )
    
    # Validate scaling if provided
    if initial_procs and hardware_type == 'gpu':
        is_valid, message = manager.validate_scaling_config(
            initial_procs, scaling_type, nodes
        )
        if not is_valid:
            logger.error(message)
            raise ValueError(f"Invalid scaling configuration:\n{message}")
        else:
            logger.info(message)
    
    # Build configuration
    config = {
        'procs_per_node': tasks_per_node,
        'gpus_per_node': gpus_per_node,
        'hardware_type': hardware_type,
    }
    
    # Add GPU-specific config
    gpu_config = manager.get_gpu_job_config()
    if gpu_config:
        config['gpu_config'] = {
            'binding': gpu_config.gpu_binding,
            'vendor': gpu_config.gpu_vendor,
            'device_var': gpu_config.device_var,
        }
    
    return config


if __name__ == '__main__':
    # Test GPU configuration
    import sys
    
    logging.basicConfig(level=logging.INFO)
    
    if len(sys.argv) < 2:
        print("Usage: python3 gpu_job_config.py <partition> [hardware_type]")
        print("Example: python3 gpu_job_config.py booster gpu")
        sys.exit(1)
    
    partition = sys.argv[1]
    hardware_type = sys.argv[2] if len(sys.argv) > 2 else 'gpu'
    
    try:
        config = configure_gpu_job(
            partition=partition,
            hardware_type=hardware_type,
            initial_procs=[2, 2, 1],
            scaling_type='weak',
            nodes=4
        )
        
        print("\n" + "="*70)
        print("GPU Job Configuration:")
        print("="*70)
        for key, value in config.items():
            print(f"{key:20s}: {value}")
        print("="*70 + "\n")
        
    except ValueError as e:
        print(f"\n❌ Configuration Error:\n{e}\n")
        sys.exit(1)
