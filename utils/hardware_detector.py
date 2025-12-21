"""
Hardware Detection Module

Automatically detects system hardware configuration:
- Cores per node
- Memory per node
- GPU count
- Architecture
"""

import os
import subprocess
import logging
from typing import Optional, Dict, Tuple

logger = logging.getLogger(__name__)


class HardwareDetector:
    """Detect hardware configuration of the system."""
    
    def __init__(self):
        self.cores_per_node = None
        self.threads_per_core = None
        self.total_memory_gb = None
        self.gpu_count = None
        self.detected = False
    
    def detect_all(self) -> Dict:
        """
        Detect all hardware parameters.
        
        Returns:
            dict: Hardware configuration
        """
        self.cores_per_node = self.detect_cores_per_node()
        self.threads_per_core = self.detect_threads_per_core()
        self.total_memory_gb = self.detect_memory()
        self.gpu_count = self.detect_gpus()
        self.detected = True
        
        return self.get_config()
    
    def detect_cores_per_node(self) -> int:
        """
        Detect number of physical CPU cores per node.
        
        Tries multiple methods:
        1. lscpu (most reliable)
        2. /proc/cpuinfo
        3. nproc
        4. SLURM environment variables
        
        Returns:
            int: Number of physical cores, or None if detection fails
        """
        # Method 1: lscpu (most reliable)
        try:
            result = subprocess.run(
                ['lscpu', '-p=core'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                cores = set()
                for line in result.stdout.strip().split('\n'):
                    if not line.startswith('#'):
                        cores.add(line.strip())
                
                core_count = len(cores)
                logger.info(f"✓ Detected {core_count} physical cores (via lscpu)")
                return core_count
        except Exception as e:
            logger.debug(f"lscpu detection failed: {e}")
        
        # Method 2: Parse /proc/cpuinfo
        try:
            with open('/proc/cpuinfo', 'r') as f:
                content = f.read()
                
            # Count unique physical core IDs
            cores = set()
            for line in content.split('\n'):
                if line.startswith('core id'):
                    core_id = line.split(':')[1].strip()
                    cores.add(core_id)
            
            if cores:
                core_count = len(cores)
                logger.info(f"✓ Detected {core_count} physical cores (via /proc/cpuinfo)")
                return core_count
        except Exception as e:
            logger.debug(f"/proc/cpuinfo detection failed: {e}")
        
        # Method 3: nproc (reports logical CPUs, may include hyperthreading)
        try:
            result = subprocess.run(
                ['nproc', '--all'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                logical_cpus = int(result.stdout.strip())
                
                # Try to detect hyperthreading
                threads_per_core = self.detect_threads_per_core()
                if threads_per_core and threads_per_core > 1:
                    core_count = logical_cpus // threads_per_core
                    logger.info(f"✓ Detected {core_count} physical cores "
                              f"({logical_cpus} logical CPUs / {threads_per_core} threads/core)")
                    return core_count
                else:
                    logger.warning(f"⚠ Detected {logical_cpus} CPUs (may include hyperthreading)")
                    return logical_cpus
        except Exception as e:
            logger.debug(f"nproc detection failed: {e}")
        
        # Method 4: SLURM environment variables
        if 'SLURM_CPUS_ON_NODE' in os.environ:
            cpus = int(os.environ['SLURM_CPUS_ON_NODE'])
            logger.info(f"✓ Detected {cpus} CPUs from SLURM_CPUS_ON_NODE")
            return cpus
        
        if 'SLURM_JOB_CPUS_PER_NODE' in os.environ:
            cpus_spec = os.environ['SLURM_JOB_CPUS_PER_NODE']
            # Parse format like "112" or "112(x2)"
            cpus = int(cpus_spec.split('(')[0])
            logger.info(f"✓ Detected {cpus} CPUs from SLURM_JOB_CPUS_PER_NODE")
            return cpus
        
        logger.warning("⚠ Could not auto-detect cores per node")
        return None
    
    def detect_threads_per_core(self) -> Optional[int]:
        """
        Detect number of hardware threads per core (hyperthreading).
        
        Returns:
            int: Threads per core (1 = no HT, 2 = HT enabled), or None
        """
        try:
            result = subprocess.run(
                ['lscpu'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                for line in result.stdout.split('\n'):
                    if 'Thread(s) per core:' in line:
                        threads = int(line.split(':')[1].strip())
                        return threads
        except Exception as e:
            logger.debug(f"Thread detection failed: {e}")
        
        return None
    
    def detect_memory(self) -> Optional[float]:
        """
        Detect total system memory in GB.
        
        Returns:
            float: Memory in GB, or None if detection fails
        """
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemTotal:'):
                        # Format: "MemTotal:       263794140 kB"
                        kb = int(line.split()[1])
                        gb = kb / (1024 * 1024)
                        logger.info(f"✓ Detected {gb:.1f} GB total memory")
                        return gb
        except Exception as e:
            logger.debug(f"Memory detection failed: {e}")
        
        return None
    
    def detect_gpus(self) -> Optional[int]:
        """
        Detect number of GPUs per node.
        
        Returns:
            int: Number of GPUs, or None if no GPUs or detection fails
        """
        # Method 1: nvidia-smi
        try:
            result = subprocess.run(
                ['nvidia-smi', '--list-gpus'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                gpu_count = len(result.stdout.strip().split('\n'))
                logger.info(f"✓ Detected {gpu_count} NVIDIA GPUs")
                return gpu_count
        except Exception as e:
            logger.debug(f"nvidia-smi detection failed: {e}")
        
        # Method 2: rocm-smi (AMD GPUs)
        try:
            result = subprocess.run(
                ['rocm-smi', '--showid'],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                lines = [l for l in result.stdout.split('\n') if 'GPU' in l]
                gpu_count = len(lines)
                if gpu_count > 0:
                    logger.info(f"✓ Detected {gpu_count} AMD GPUs")
                    return gpu_count
        except Exception as e:
            logger.debug(f"rocm-smi detection failed: {e}")
        
        # Method 3: Check /dev
        try:
            import glob
            nvidia_devices = glob.glob('/dev/nvidia[0-9]*')
            if nvidia_devices:
                gpu_count = len(nvidia_devices)
                logger.info(f"✓ Detected {gpu_count} GPU devices in /dev")
                return gpu_count
        except Exception as e:
            logger.debug(f"/dev detection failed: {e}")
        
        logger.debug("No GPUs detected")
        return 0
    
    def get_config(self) -> Dict:
        """Get detected hardware configuration."""
        return {
            'cores_per_node': self.cores_per_node,
            'threads_per_core': self.threads_per_core,
            'total_memory_gb': self.total_memory_gb,
            'gpu_count': self.gpu_count,
            'detected': self.detected
        }
    
    def print_summary(self):
        """Print hardware detection summary."""
        print()
        print("="*70)
        print("Hardware Detection Summary")
        print("="*70)
        
        if self.cores_per_node:
            print(f"  Cores per node:    {self.cores_per_node}")
        else:
            print(f"  Cores per node:    NOT DETECTED")
        
        if self.threads_per_core:
            print(f"  Threads per core:  {self.threads_per_core}")
            if self.cores_per_node and self.threads_per_core:
                logical = self.cores_per_node * self.threads_per_core
                print(f"  Logical CPUs:      {logical}")
        
        if self.total_memory_gb:
            print(f"  Total memory:      {self.total_memory_gb:.1f} GB")
        
        if self.gpu_count is not None:
            if self.gpu_count > 0:
                print(f"  GPUs per node:     {self.gpu_count}")
            else:
                print(f"  GPUs per node:     0 (CPU-only node)")
        
        print("="*70)
        print()


def detect_hardware() -> Dict:
    """
    Convenience function to detect hardware and return config.
    
    Returns:
        dict: Hardware configuration
    """
    detector = HardwareDetector()
    config = detector.detect_all()
    detector.print_summary()
    return config


if __name__ == "__main__":
    # Test hardware detection
    logging.basicConfig(level=logging.INFO)
    config = detect_hardware()
    
    print("Detected configuration:")
    for key, value in config.items():
        print(f"  {key}: {value}")
