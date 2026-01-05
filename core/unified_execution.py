#!/usr/bin/env python3
"""
Unified Execution Module for HPC-ScaleTest

This module provides the SINGLE SOURCE OF TRUTH for:
1. Hardware topology detection (CPUs/GPUs per node)
2. MPI mapping computation (ranks per node, cores per rank)
3. MPI command generation (mpirun -np X --map-by ppr:Y:node:PE=Z)
4. GPU binding script generation (bind.sh)

Design Principles:
==================
- NO HARDCODED VALUES - All topology is detected at runtime
- SYSTEM AGNOSTIC - Works on Leonardo (NVIDIA), LUMI (AMD), and other systems
- GPU VENDOR DETECTION - Auto-detects NVIDIA via nvidia-smi, AMD via rocm-smi
- SLURM INTEGRATION - Uses SLURM environment variables when available

Usage:
======
    executor = UnifiedExecutor(partition="boost_usr_prod")
    
    # Auto-detect topology
    topology = executor.detect_topology()
    
    # Generate MPI command for N nodes
    mpi_cmd = executor.generate_mpi_command(
        num_nodes=4,
        executable="$BINARY/iPIC3D",
        args=["os-stdin"]
    )
    
    # Generate complete job script
    script = executor.generate_job_script(
        num_nodes=4,
        executable="$BINARY/iPIC3D",
        args=["os-stdin"],
        job_name="ipic3d_test"
    )

Author: HPC-ScaleTest Contributors
"""

import os
import subprocess
import logging
import re
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Any
from enum import Enum
from pathlib import Path

logger = logging.getLogger(__name__)


class GPUVendor(Enum):
    """Supported GPU vendors."""
    NVIDIA = "nvidia"
    AMD = "amd"
    INTEL = "intel"
    NONE = "none"


@dataclass
class HardwareTopology:
    """
    Complete hardware topology for a compute node.
    
    All values are detected at runtime - NO hardcoded defaults.
    
    Attributes:
        cpu_cores_per_node: Total CPU cores available per node
        gpus_per_node: Number of GPUs per node (0 for CPU-only)
        gpu_vendor: Detected GPU vendor (nvidia, amd, intel, none)
        gpu_model: GPU model name if detected
        partition: SLURM partition name
        detection_method: How the topology was detected
        
    Derived (computed in __post_init__):
        ranks_per_node: MPI ranks per node (= GPUs for GPU jobs, = cores for CPU jobs)
        cores_per_rank: CPU cores allocated per MPI rank
    """
    cpu_cores_per_node: int
    gpus_per_node: int = 0
    gpu_vendor: GPUVendor = GPUVendor.NONE
    gpu_model: str = ""
    partition: str = ""
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
            # CPU job: 1 MPI rank per core
            self.ranks_per_node = self.cpu_cores_per_node
            self.cores_per_rank = 1
    
    def __str__(self) -> str:
        return (
            f"HardwareTopology(partition={self.partition}, "
            f"cpus={self.cpu_cores_per_node}, gpus={self.gpus_per_node}, "
            f"vendor={self.gpu_vendor.value}, "
            f"ranks/node={self.ranks_per_node}, cores/rank={self.cores_per_rank}, "
            f"method={self.detection_method})"
        )


class TopologyDetector:
    """
    Unified hardware topology detector.
    
    Detection Strategy (priority order):
    1. SLURM environment variables (inside a job)
    2. sinfo/scontrol queries (on login node)
    3. Direct system introspection (nvidia-smi, rocm-smi, lscpu)
    
    GPU Detection:
    - NVIDIA: nvidia-smi --query-gpu=count --format=csv,noheader
    - AMD: rocm-smi --showproductname (count GPU entries)
    """
    
    def __init__(self):
        self._cache: Dict[str, HardwareTopology] = {}
    
    def detect(self, partition: Optional[str] = None) -> HardwareTopology:
        """
        Detect hardware topology for a partition.
        
        Args:
            partition: SLURM partition name. If None, detects from environment.
        
        Returns:
            HardwareTopology with all detected values
        
        Raises:
            RuntimeError: If topology cannot be detected
        """
        # Resolve partition from environment if not provided
        if not partition:
            partition = os.environ.get('SLURM_JOB_PARTITION', '')
        
        # Check cache
        cache_key = partition or '_default_'
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # Detection chain
        topology = None
        method = ""
        
        # Method 1: SLURM environment (highest priority if inside a job)
        if os.environ.get('SLURM_JOB_ID'):
            topology = self._detect_from_slurm_env(partition)
            if topology:
                method = "SLURM environment variables"
        
        # Method 2: sinfo query (on login node)
        if topology is None and partition:
            topology = self._detect_from_sinfo(partition)
            if topology:
                method = "sinfo query"
        
        # Method 3: scontrol query
        if topology is None and partition:
            topology = self._detect_from_scontrol(partition)
            if topology:
                method = "scontrol query"
        
        # Method 4: Direct system introspection
        if topology is None:
            topology = self._detect_from_system()
            if topology:
                method = "system introspection"
        
        if topology is None:
            raise RuntimeError(
                f"Cannot detect topology for partition '{partition}'. "
                "Ensure you're on an HPC system with SLURM access."
            )
        
        topology.detection_method = method
        topology.partition = partition or ""
        
        # Cache and return
        self._cache[cache_key] = topology
        
        logger.info(f"Detected topology: {topology}")
        return topology
    
    def _detect_from_slurm_env(self, partition: str) -> Optional[HardwareTopology]:
        """Detect topology from SLURM environment variables."""
        try:
            cpu_cores = None
            gpus = 0
            gpu_vendor = GPUVendor.NONE
            
            # CPU cores
            if 'SLURM_CPUS_ON_NODE' in os.environ:
                cpu_cores = int(os.environ['SLURM_CPUS_ON_NODE'])
            elif 'SLURM_JOB_CPUS_PER_NODE' in os.environ:
                # Format: "32" or "32(x2)" for heterogeneous
                val = os.environ['SLURM_JOB_CPUS_PER_NODE']
                cpu_cores = int(val.split('(')[0])
            
            if cpu_cores is None:
                return None
            
            # GPUs - try multiple sources
            if 'SLURM_GPUS_ON_NODE' in os.environ:
                gpus = int(os.environ['SLURM_GPUS_ON_NODE'])
            elif 'SLURM_GPUS_PER_NODE' in os.environ:
                gpus = int(os.environ['SLURM_GPUS_PER_NODE'])
            elif 'CUDA_VISIBLE_DEVICES' in os.environ:
                devs = os.environ['CUDA_VISIBLE_DEVICES']
                if devs and devs.lower() not in ('', 'nodevfiles'):
                    gpus = len([d for d in devs.split(',') if d.strip()])
                    gpu_vendor = GPUVendor.NVIDIA
            elif 'ROCR_VISIBLE_DEVICES' in os.environ:
                devs = os.environ['ROCR_VISIBLE_DEVICES']
                if devs:
                    gpus = len([d for d in devs.split(',') if d.strip()])
                    gpu_vendor = GPUVendor.AMD
            
            # Detect GPU vendor if GPUs found but vendor unknown
            if gpus > 0 and gpu_vendor == GPUVendor.NONE:
                gpu_vendor = self._detect_gpu_vendor()
            
            return HardwareTopology(
                cpu_cores_per_node=cpu_cores,
                gpus_per_node=gpus,
                gpu_vendor=gpu_vendor,
                partition=partition
            )
        except Exception as e:
            logger.debug(f"SLURM env detection failed: {e}")
            return None
    
    def _detect_from_sinfo(self, partition: str) -> Optional[HardwareTopology]:
        """Detect topology using sinfo command."""
        try:
            # Query CPUs and GRES (GPU resources)
            result = subprocess.run(
                ['sinfo', '-p', partition, '-N', '-h', '-o', '%c %G'],
                capture_output=True, text=True, timeout=10
            )
            
            if result.returncode != 0 or not result.stdout.strip():
                return None
            
            # Parse first line
            line = result.stdout.strip().split('\n')[0]
            parts = line.split()
            
            cpu_cores = None
            gpus = 0
            gpu_vendor = GPUVendor.NONE
            
            for part in parts:
                # CPU count is usually first numeric value
                if cpu_cores is None and part.isdigit():
                    cpu_cores = int(part)
                # GRES format: gpu:model:count or gpu:count
                elif 'gpu' in part.lower():
                    match = re.search(r'gpu:(?:(\w+):)?(\d+)', part, re.IGNORECASE)
                    if match:
                        gpu_model = match.group(1) or ""
                        gpus = int(match.group(2))
                        
                        # Infer vendor from model name
                        if gpu_model:
                            gpu_model_lower = gpu_model.lower()
                            if any(n in gpu_model_lower for n in ['a100', 'v100', 'h100', 'a40', 'rtx', 'tesla']):
                                gpu_vendor = GPUVendor.NVIDIA
                            elif any(a in gpu_model_lower for a in ['mi100', 'mi250', 'mi300', 'instinct']):
                                gpu_vendor = GPUVendor.AMD
            
            if cpu_cores is None:
                return None
            
            return HardwareTopology(
                cpu_cores_per_node=cpu_cores,
                gpus_per_node=gpus,
                gpu_vendor=gpu_vendor,
                partition=partition
            )
        except Exception as e:
            logger.debug(f"sinfo detection failed: {e}")
            return None
    
    def _detect_from_scontrol(self, partition: str) -> Optional[HardwareTopology]:
        """Detect topology using scontrol command."""
        try:
            result = subprocess.run(
                ['scontrol', 'show', 'partition', partition],
                capture_output=True, text=True, timeout=10
            )
            
            if result.returncode != 0:
                return None
            
            output = result.stdout
            cpu_cores = None
            gpus = 0
            gpu_vendor = GPUVendor.NONE
            
            # Parse MaxCPUsPerNode
            match = re.search(r'MaxCPUsPerNode=(\d+)', output)
            if match:
                cpu_cores = int(match.group(1))
            
            # Parse GRES
            match = re.search(r'GRES=([^\s]+)', output)
            if match:
                gres = match.group(1)
                gpu_match = re.search(r'gpu:(?:(\w+):)?(\d+)', gres)
                if gpu_match:
                    gpu_model = gpu_match.group(1) or ""
                    gpus = int(gpu_match.group(2))
            
            # If no CPU count from partition, try querying a node
            if cpu_cores is None:
                match = re.search(r'Nodes=([^\s,\[]+)', output)
                if match:
                    node_result = subprocess.run(
                        ['scontrol', 'show', 'node', match.group(1)],
                        capture_output=True, text=True, timeout=10
                    )
                    if node_result.returncode == 0:
                        cpu_match = re.search(r'CPUTot=(\d+)', node_result.stdout)
                        if cpu_match:
                            cpu_cores = int(cpu_match.group(1))
            
            if cpu_cores is None:
                return None
            
            return HardwareTopology(
                cpu_cores_per_node=cpu_cores,
                gpus_per_node=gpus,
                gpu_vendor=gpu_vendor,
                partition=partition
            )
        except Exception as e:
            logger.debug(f"scontrol detection failed: {e}")
            return None
    
    def _detect_from_system(self) -> Optional[HardwareTopology]:
        """Detect topology from direct system introspection."""
        try:
            # CPU cores via lscpu or os.cpu_count()
            cpu_cores = self._detect_cpu_cores()
            if cpu_cores is None:
                return None
            
            # GPU detection - try NVIDIA first, then AMD
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
            
            return HardwareTopology(
                cpu_cores_per_node=cpu_cores,
                gpus_per_node=gpus,
                gpu_vendor=gpu_vendor,
                gpu_model=gpu_model
            )
        except Exception as e:
            logger.debug(f"System introspection failed: {e}")
            return None
    
    def _detect_cpu_cores(self) -> Optional[int]:
        """Detect CPU cores using lscpu or os.cpu_count()."""
        try:
            result = subprocess.run(
                ['lscpu', '--parse=CPU'],
                capture_output=True, text=True, timeout=5
            )
            if result.returncode == 0:
                # Count non-comment lines
                lines = [l for l in result.stdout.strip().split('\n') if not l.startswith('#')]
                return len(lines)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        # Fallback to Python
        return os.cpu_count()
    
    def _detect_nvidia_gpus(self) -> Optional[Dict]:
        """Detect NVIDIA GPUs using nvidia-smi."""
        try:
            result = subprocess.run(
                ['nvidia-smi', '--query-gpu=count,name', '--format=csv,noheader,nounits'],
                capture_output=True, text=True, timeout=10
            )
            
            if result.returncode == 0 and result.stdout.strip():
                lines = result.stdout.strip().split('\n')
                count = len(lines)
                
                # Parse model from first line
                model = ""
                if lines:
                    parts = lines[0].split(',')
                    if len(parts) > 1:
                        model = parts[1].strip()
                
                return {'count': count, 'model': model}
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        # Fallback: check /dev/nvidia*
        try:
            import glob
            nvidia_devs = glob.glob('/dev/nvidia[0-9]*')
            if nvidia_devs:
                return {'count': len(nvidia_devs), 'model': ''}
        except Exception:
            pass
        
        return None
    
    def _detect_amd_gpus(self) -> Optional[Dict]:
        """Detect AMD GPUs using rocm-smi."""
        try:
            result = subprocess.run(
                ['rocm-smi', '--showproductname'],
                capture_output=True, text=True, timeout=10
            )
            
            if result.returncode == 0:
                # Count GPU entries
                count = len(re.findall(r'GPU\[(\d+)\]', result.stdout))
                if count > 0:
                    model_match = re.search(r'Card series:\s*(.+)', result.stdout)
                    model = model_match.group(1).strip() if model_match else ''
                    return {'count': count, 'model': model}
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        return None
    
    def _detect_gpu_vendor(self) -> GPUVendor:
        """Detect GPU vendor when GPUs are present but vendor unknown."""
        if self._detect_nvidia_gpus():
            return GPUVendor.NVIDIA
        if self._detect_amd_gpus():
            return GPUVendor.AMD
        return GPUVendor.NONE


class UnifiedExecutor:
    """
    Unified execution engine for HPC jobs.
    
    This class provides:
    1. Automatic hardware topology detection
    2. Correct MPI command generation
    3. GPU binding script generation
    4. Complete job script generation
    
    Example:
        executor = UnifiedExecutor(partition="boost_usr_prod")
        
        # For Leonardo Booster (32 cores, 4 GPUs per node):
        # - 4 MPI ranks per node (1 per GPU)
        # - 8 CPU cores per rank (32/4)
        
        # 1 node: mpirun -np 4 --map-by ppr:4:node:PE=8 ./bind.sh app
        # 2 nodes: mpirun -np 8 --map-by ppr:4:node:PE=8 ./bind.sh app
        # 4 nodes: mpirun -np 16 --map-by ppr:4:node:PE=8 ./bind.sh app
    """
    
    def __init__(
        self,
        partition: str = "",
        account: str = "",
        qos: str = "",
        time_limit: str = "02:00:00"
    ):
        """
        Initialize the unified executor.
        
        Args:
            partition: SLURM partition name
            account: SLURM account
            qos: Quality of service
            time_limit: Default time limit for jobs
        """
        self.partition = partition
        self.account = account
        self.qos = qos
        self.time_limit = time_limit
        
        self._detector = TopologyDetector()
        self._topology: Optional[HardwareTopology] = None
    
    def detect_topology(self) -> HardwareTopology:
        """
        Detect hardware topology for the partition.
        
        Returns:
            HardwareTopology with all detected values
        """
        if self._topology is None:
            self._topology = self._detector.detect(self.partition)
        return self._topology
    
    @property
    def topology(self) -> HardwareTopology:
        """Get topology (auto-detect if needed)."""
        return self.detect_topology()
    
    def generate_mpi_command(
        self,
        num_nodes: int,
        executable: str,
        args: List[str] = None,
        bind_script: str = "./bind.sh",
        report_bindings: bool = True
    ) -> str:
        """
        Generate mpirun command with correct topology-based arguments.
        
        Algorithm:
        - total_ranks = num_nodes * ranks_per_node
        - For GPU jobs: ranks_per_node = GPUs per node
        - cores_per_rank = CPU cores per node / ranks_per_node
        - Command: mpirun -np <total> --map-by ppr:<ranks>:node:PE=<cores> [--report-bindings] [bind.sh] <exe> [args]
        
        Args:
            num_nodes: Number of nodes
            executable: Path to executable (e.g., "$BINARY/iPIC3D")
            args: Arguments for executable
            bind_script: Path to GPU binding script
            report_bindings: Include --report-bindings flag
        
        Returns:
            Complete mpirun command string
        
        Example (Leonardo Booster, 4 nodes):
            mpirun -np 16 --map-by ppr:4:node:PE=8 --report-bindings ./bind.sh $BINARY/iPIC3D os-stdin
        """
        args = args or []
        topo = self.topology
        
        # Calculate MPI layout
        total_ranks = num_nodes * topo.ranks_per_node
        
        # Build command parts
        parts = [
            'mpirun',
            '-np', str(total_ranks),
            '--map-by', f'ppr:{topo.ranks_per_node}:node:PE={topo.cores_per_rank}'
        ]
        
        if report_bindings:
            parts.append('--report-bindings')
        
        # GPU binding wrapper (for GPU jobs)
        if topo.gpus_per_node > 0 and bind_script:
            parts.append(bind_script)
        
        # Executable and arguments
        parts.append(executable)
        parts.extend(args)
        
        return ' '.join(parts)
    
    def generate_gpu_binding_script(self) -> str:
        """
        Generate GPU binding script content.
        
        The script sets CUDA_VISIBLE_DEVICES (NVIDIA) or ROCR_VISIBLE_DEVICES (AMD)
        based on the MPI local rank (OMPI_COMM_WORLD_LOCAL_RANK).
        
        Returns:
            Shell script content as string
        """
        topo = self.topology
        
        # Select environment variable based on GPU vendor
        if topo.gpu_vendor == GPUVendor.AMD:
            gpu_env_var = "ROCR_VISIBLE_DEVICES"
        elif topo.gpu_vendor == GPUVendor.INTEL:
            gpu_env_var = "ZE_AFFINITY_MASK"
        else:
            # Default to NVIDIA
            gpu_env_var = "CUDA_VISIBLE_DEVICES"
        
        return f'''#!/bin/bash
# GPU Binding Script - Generated by HPC-ScaleTest
# 
# This script binds each MPI rank to a unique GPU based on local rank.
# GPU Vendor: {topo.gpu_vendor.value}
# GPUs per node: {topo.gpus_per_node}
# Ranks per node: {topo.ranks_per_node}

# Determine local rank from available MPI environment variables
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

# Bind this rank to its assigned GPU
export {gpu_env_var}=$LOCAL_RANK

# Execute the application
exec "$@"
'''
    
    def generate_slurm_directives(
        self,
        num_nodes: int,
        job_name: str,
        time_limit: str = None,
        account: str = None,
        qos: str = None,
        output_dir: str = ".",
        exclusive: bool = True
    ) -> List[str]:
        """
        Generate SLURM directives for job script.
        
        IMPORTANT: Uses FULL node resources (all CPUs) for --ntasks-per-node
        to ensure proper resource allocation. MPI execution uses different
        (usually smaller) task count.
        
        Args:
            num_nodes: Number of nodes
            job_name: Job name
            time_limit: Override time limit
            account: Override account
            qos: Override QoS
            output_dir: Output directory for logs
            exclusive: Request exclusive node access
        
        Returns:
            List of SBATCH directive lines
        """
        topo = self.topology
        time_limit = time_limit or self.time_limit
        account = account or self.account
        qos = qos or self.qos
        
        directives = [
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --nodes={num_nodes}",
            # CRITICAL: Use FULL CPU count for resource allocation
            f"#SBATCH --ntasks-per-node={topo.cpu_cores_per_node}",
        ]
        
        if self.partition:
            directives.append(f"#SBATCH --partition={self.partition}")
        
        if topo.gpus_per_node > 0:
            directives.append(f"#SBATCH --gres=gpu:{topo.gpus_per_node}")
        
        if account:
            directives.append(f"#SBATCH --account={account}")
        
        if qos:
            directives.append(f"#SBATCH --qos={qos}")
        
        directives.append(f"#SBATCH --time={time_limit}")
        directives.append(f"#SBATCH --output={output_dir}/job_%j.out")
        directives.append(f"#SBATCH --error={output_dir}/job_%j.err")
        
        if exclusive:
            directives.append("#SBATCH --exclusive")
        
        return directives
    
    def generate_job_script(
        self,
        num_nodes: int,
        executable: str,
        args: List[str] = None,
        job_name: str = "hpc_job",
        time_limit: str = None,
        account: str = None,
        qos: str = None,
        output_dir: str = ".",
        env_setup: List[str] = None,
        binary_dir: str = "",
        report_bindings: bool = True
    ) -> str:
        """
        Generate complete SLURM job script.
        
        Args:
            num_nodes: Number of nodes
            executable: Executable name (e.g., "iPIC3D")
            args: Arguments for executable
            job_name: Job name
            time_limit: Time limit
            account: SLURM account
            qos: Quality of service
            output_dir: Working directory and output location
            env_setup: Environment setup commands (module loads, exports)
            binary_dir: Directory containing the executable
            report_bindings: Include --report-bindings in mpirun
        
        Returns:
            Complete job script as string
        """
        args = args or []
        env_setup = env_setup or []
        topo = self.topology
        
        # Build executable path
        if binary_dir:
            exe_path = f"$BINARY/{executable}"
        else:
            exe_path = executable
        
        # Generate mpirun command
        mpi_cmd = self.generate_mpi_command(
            num_nodes=num_nodes,
            executable=exe_path,
            args=args,
            bind_script="./bind.sh" if topo.gpus_per_node > 0 else None,
            report_bindings=report_bindings
        )
        
        # Generate SLURM directives
        directives = self.generate_slurm_directives(
            num_nodes=num_nodes,
            job_name=job_name,
            time_limit=time_limit,
            account=account,
            qos=qos,
            output_dir=output_dir
        )
        
        # Generate GPU binding script content
        bind_script_content = self.generate_gpu_binding_script()
        
        # Calculate total MPI ranks for documentation
        total_ranks = num_nodes * topo.ranks_per_node
        
        # Build script
        script = f'''#!/bin/bash
# =============================================================================
# SLURM Job Script - Generated by HPC-ScaleTest
# =============================================================================
# Partition: {self.partition}
# Detection Method: {topo.detection_method}
#
# Hardware Topology (per node):
#   CPU cores: {topo.cpu_cores_per_node}
#   GPUs: {topo.gpus_per_node}
#   GPU vendor: {topo.gpu_vendor.value}
#
# MPI Mapping (computed from topology):
#   Ranks per node: {topo.ranks_per_node} ({"1 per GPU" if topo.gpus_per_node > 0 else "1 per core"})
#   Cores per rank: {topo.cores_per_rank}
#   Total ranks: {total_ranks}
#
# IMPORTANT: SLURM resource allocation vs MPI execution
# -----------------------------------------------------
# SLURM allocates FULL node resources: --ntasks-per-node={topo.cpu_cores_per_node}
# MPI executes with derived tasks: -np {total_ranks}
# =============================================================================

{chr(10).join(directives)}

# =============================================================================
# Environment Setup
# =============================================================================
'''
        
        for cmd in env_setup:
            script += f"{cmd}\n"
        
        script += f'''
# Binary location
export BINARY={binary_dir if binary_dir else output_dir}

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_NNODES"
echo "Partition: $SLURM_JOB_PARTITION"
echo "CPUs on node: $SLURM_CPUS_ON_NODE"
echo "GPUs on node: ${{SLURM_GPUS_ON_NODE:-{topo.gpus_per_node}}}"
echo "BINARY: $BINARY"
echo "=========================================="

cd {output_dir}

# Set OpenMP threads to cores per rank
export OMP_NUM_THREADS={topo.cores_per_rank}
'''
        
        # Add GPU binding script if GPU job
        if topo.gpus_per_node > 0:
            script += f'''
# =============================================================================
# GPU Binding Script
# =============================================================================
cat > bind.sh << 'BIND_EOF'
{bind_script_content}BIND_EOF
chmod +x bind.sh
'''
        
        script += f'''
# =============================================================================
# Execute Application
# =============================================================================
echo ""
echo "Starting: $(date)"
echo "Command: {mpi_cmd}"
echo ""

start_time=$(date +%s.%N)

{mpi_cmd}

exit_code=$?
end_time=$(date +%s.%N)
elapsed=$(echo "$end_time - $start_time" | bc)

echo ""
echo "=========================================="
echo "Completed: $(date)"
echo "Exit code: $exit_code"
echo "Elapsed: $elapsed seconds"
echo "=========================================="

exit $exit_code
'''
        
        return script


# =============================================================================
# Global Singleton
# =============================================================================

_global_executor: Optional[UnifiedExecutor] = None


def get_executor(partition: str = "") -> UnifiedExecutor:
    """Get or create global executor instance."""
    global _global_executor
    if _global_executor is None or (partition and _global_executor.partition != partition):
        _global_executor = UnifiedExecutor(partition=partition)
    return _global_executor


def detect_topology(partition: str = "") -> HardwareTopology:
    """Convenience function to detect topology."""
    return get_executor(partition).detect_topology()


def generate_mpi_command(
    partition: str,
    num_nodes: int,
    executable: str,
    args: List[str] = None
) -> str:
    """Convenience function to generate MPI command."""
    executor = get_executor(partition)
    return executor.generate_mpi_command(num_nodes, executable, args)


# =============================================================================
# Self-Test
# =============================================================================

if __name__ == '__main__':
    import sys
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s - %(message)s'
    )
    
    print("=" * 70)
    print(" Unified Execution Module - Self Test")
    print("=" * 70)
    
    # Test with mock topology (Leonardo Booster: 32 cores, 4 GPUs)
    print("\n[Test 1] Mock Leonardo Booster (32 cores, 4 GPUs)")
    
    mock_topo = HardwareTopology(
        cpu_cores_per_node=32,
        gpus_per_node=4,
        gpu_vendor=GPUVendor.NVIDIA,
        partition="boost_usr_prod",
        detection_method="mock"
    )
    
    print(f"  Topology: {mock_topo}")
    print(f"  Ranks per node: {mock_topo.ranks_per_node} (expected: 4)")
    print(f"  Cores per rank: {mock_topo.cores_per_rank} (expected: 8)")
    
    assert mock_topo.ranks_per_node == 4
    assert mock_topo.cores_per_rank == 8
    
    # Create executor with mock topology
    executor = UnifiedExecutor(partition="boost_usr_prod")
    executor._topology = mock_topo
    
    print("\n[Test 2] MPI command generation")
    for nodes in [1, 2, 4]:
        cmd = executor.generate_mpi_command(
            num_nodes=nodes,
            executable="$BINARY/iPIC3D",
            args=["os-stdin"]
        )
        expected_np = nodes * 4
        print(f"  {nodes} node(s): {cmd}")
        assert f"-np {expected_np}" in cmd
        assert "ppr:4:node:PE=8" in cmd
        assert "./bind.sh" in cmd
    
    print("\n[Test 3] GPU binding script")
    bind_script = executor.generate_gpu_binding_script()
    print("  Generated bind.sh content:")
    for line in bind_script.split('\n')[:10]:
        print(f"    {line}")
    assert "CUDA_VISIBLE_DEVICES" in bind_script
    assert "OMPI_COMM_WORLD_LOCAL_RANK" in bind_script
    
    print("\n[Test 4] SLURM directives")
    directives = executor.generate_slurm_directives(
        num_nodes=4,
        job_name="test_job",
        account="my_account"
    )
    directives_str = "\n".join(directives)
    print("  Generated directives:")
    for d in directives[:5]:
        print(f"    {d}")
    
    assert "--ntasks-per-node=32" in directives_str, "Should use FULL CPU count"
    assert "--gres=gpu:4" in directives_str, "Should request GPUs"
    
    print("\n[Test 5] Complete job script")
    script = executor.generate_job_script(
        num_nodes=4,
        executable="iPIC3D",
        args=["os-stdin"],
        job_name="ipic3d_test",
        binary_dir="/path/to/build",
        env_setup=["module load cuda/12.0"]
    )
    
    # Verify key elements
    checks = [
        ("mpirun -np 16" in script, "mpirun -np 16 (4 nodes × 4 GPUs)"),
        ("ppr:4:node:PE=8" in script, "--map-by ppr:4:node:PE=8"),
        ("./bind.sh $BINARY/iPIC3D" in script, "./bind.sh wrapping executable"),
        ("CUDA_VISIBLE_DEVICES=$LOCAL_RANK" in script, "GPU binding in bind.sh"),
        ("--ntasks-per-node=32" in script, "SLURM full node allocation"),
    ]
    
    print("  Verification:")
    all_passed = True
    for passed, desc in checks:
        status = "✓" if passed else "✗"
        print(f"    {status} {desc}")
        if not passed:
            all_passed = False
    
    # Test AMD GPU support
    print("\n[Test 6] AMD GPU support (LUMI-style: 8 GPUs)")
    amd_topo = HardwareTopology(
        cpu_cores_per_node=64,
        gpus_per_node=8,
        gpu_vendor=GPUVendor.AMD,
        partition="gpu",
        detection_method="mock"
    )
    
    amd_executor = UnifiedExecutor(partition="gpu")
    amd_executor._topology = amd_topo
    
    bind_script = amd_executor.generate_gpu_binding_script()
    print(f"  GPU env var: {'ROCR_VISIBLE_DEVICES' if 'ROCR_VISIBLE_DEVICES' in bind_script else 'UNKNOWN'}")
    assert "ROCR_VISIBLE_DEVICES" in bind_script, "Should use AMD GPU env var"
    
    cmd = amd_executor.generate_mpi_command(2, "app", ["input"])
    print(f"  2 nodes: {cmd}")
    assert "-np 16" in cmd  # 2 nodes × 8 GPUs
    assert "ppr:8:node:PE=8" in cmd  # 8 ranks, 64/8=8 cores each
    
    print("\n" + "=" * 70)
    if all_passed:
        print(" All tests PASSED ✓")
    else:
        print(" Some tests FAILED ✗")
        sys.exit(1)
    print("=" * 70)
