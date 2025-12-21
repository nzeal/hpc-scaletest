#!/usr/bin/env python3
"""
Advanced GPU Resource Manager - Python 3.6+ Compatible

Handles the distinction between:
- CPU allocation (--ntasks-per-node for SLURM)
- GPU allocation (--gres=gpu:N)
- Actual MPI tasks (1 per GPU)
- Cores per task (total_cores / num_gpus)

Compatible with Python 3.6+ (no dataclasses, no type hints in function signatures)
MPI-agnostic (works with OpenMPI, Intel MPI, MPICH, etc.)
"""

import logging
import subprocess
import re
import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from utils.mpi_detector import MPIDetector, MPIInfo
    HAS_MPI_DETECTOR = True
except ImportError:
    HAS_MPI_DETECTOR = False
    logger = logging.getLogger(__name__)
    logger.warning("MPI detector not available, using safe defaults")


logger = logging.getLogger(__name__)


class GPUNodeConfig(object):
    """Complete GPU node configuration."""
    
    def __init__(self, gpus_per_node, cpus_per_node, cores_per_gpu,
                 memory_per_node_gb, gpu_model, gpu_vendor, partition):
        self.gpus_per_node = gpus_per_node
        self.cpus_per_node = cpus_per_node
        self.cores_per_gpu = cores_per_gpu
        self.memory_per_node_gb = memory_per_node_gb
        self.gpu_model = gpu_model
        self.gpu_vendor = gpu_vendor
        self.partition = partition
    
    def __repr__(self):
        return (
            "GPUNodeConfig(gpus_per_node={}, cpus_per_node={}, "
            "cores_per_gpu={}, gpu_model='{}', partition='{}')"
        ).format(
            self.gpus_per_node, self.cpus_per_node, self.cores_per_gpu,
            self.gpu_model, self.partition
        )


class AdvancedGPUManager(object):
    """
    Manages GPU resources with full awareness of CPU/GPU topology.
    
    Distinguishes between:
    - CPU cores available (e.g., 32 on Leonardo Booster)
    - GPU count (e.g., 4 on Leonardo Booster)
    - MPI tasks = GPU count (1:1 mapping)
    - Cores per MPI task = total_cores / gpu_count
    
    Supports flexible partition names (both 'booster' and 'boost_usr_prod')
    """
    
    # Known partition name mappings for different systems
    PARTITION_ALIASES = {
        'booster': ['booster', 'boost_usr_prod', 'boost_qos_dbg', 'boost_qos_bprod'],
        'dcgp': ['dcgp', 'dcgp_usr_prod', 'dcgp_usr_preempt'],
    }
    
    def __init__(self):
        self.node_config = None
    
    def _normalize_partition(self, partition):
        """
        Normalize partition name to handle aliases.
        
        Examples:
        - 'boost_usr_prod' -> 'booster'
        - 'booster' -> 'booster'
        - 'boost_qos_dbg' -> 'booster'
        """
        partition_lower = partition.lower()
        
        # Check if it matches any known aliases
        for canonical, aliases in self.PARTITION_ALIASES.items():
            if partition_lower in [a.lower() for a in aliases]:
                logger.debug("Normalized partition '{}' -> '{}'".format(partition, canonical))
                return canonical
        
        return partition
    
    def _find_partition(self, partition_hint):
        """
        Find actual partition name from hint.
        
        Tries:
        1. Exact match
        2. Alias match
        3. Fuzzy match (contains hint)
        """
        try:
            # Get list of available partitions
            cmd = ["sinfo", "-o", "%R", "-h"]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode != 0:
                return partition_hint
            
            available_partitions = [p.strip() for p in result.stdout.split('\n') if p.strip()]
            
            # Try exact match first
            if partition_hint in available_partitions:
                logger.debug("Found exact partition match: '{}'".format(partition_hint))
                return partition_hint
            
            # Try alias matching
            normalized = self._normalize_partition(partition_hint)
            for canonical, aliases in self.PARTITION_ALIASES.items():
                if normalized == canonical:
                    # Find which alias exists on this system
                    for alias in aliases:
                        if alias in available_partitions:
                            logger.info("✓ Mapped partition '{}' -> '{}' (actual name on system)".format(
                                partition_hint, alias))
                            return alias
            
            # Try fuzzy match (contains)
            partition_lower = partition_hint.lower()
            for part in available_partitions:
                if partition_lower in part.lower() or part.lower() in partition_lower:
                    logger.info("✓ Fuzzy matched partition '{}' -> '{}'".format(partition_hint, part))
                    return part
            
            # No match, return original
            logger.warning("⚠ Could not find partition '{}' on system".format(partition_hint))
            logger.warning("   Available partitions: {}".format(', '.join(available_partitions[:5])))
            return partition_hint
            
        except Exception as e:
            logger.debug("Partition lookup failed: {}".format(e))
            return partition_hint
    
    def detect_gpu_node_config(self, partition):
        """
        Detect complete GPU node configuration from SLURM partition.
        
        Supports flexible partition naming (e.g., 'booster' or 'boost_usr_prod').
        """
        # Find actual partition name
        actual_partition = self._find_partition(partition)
        
        try:
            # Query partition information
            cmd = ["scontrol", "show", "partition", actual_partition]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode != 0:
                logger.warning("Could not query partition '{}'".format(actual_partition))
                return None
            
            partition_info = result.stdout
            
            # Parse GPU count
            gpus_per_node = self._parse_gres_gpu(partition_info)
            if not gpus_per_node:
                logger.warning("No GPU information found for partition '{}'".format(actual_partition))
                return None
            
            # Query node information from a sample node
            node_name = self._get_sample_node(actual_partition)
            if node_name:
                cpus_per_node, memory_gb = self._query_node_resources(node_name)
            else:
                # Fallback to partition-level info
                cpus_per_node = self._parse_cpus_per_node(partition_info)
                memory_gb = 0
            
            # Calculate cores per GPU
            cores_per_gpu = cpus_per_node // gpus_per_node if gpus_per_node > 0 else 1
            
            # Detect GPU model
            gpu_model, gpu_vendor = self._detect_gpu_model(node_name)
            
            config = GPUNodeConfig(
                gpus_per_node=gpus_per_node,
                cpus_per_node=cpus_per_node,
                cores_per_gpu=cores_per_gpu,
                memory_per_node_gb=memory_gb,
                gpu_model=gpu_model,
                gpu_vendor=gpu_vendor,
                partition=actual_partition
            )
            
            self.node_config = config
            
            logger.info("✓ Detected GPU node configuration:")
            logger.info("  Partition: {}".format(actual_partition))
            logger.info("  GPUs per node: {} ({})".format(gpus_per_node, gpu_model))
            logger.info("  CPUs per node: {}".format(cpus_per_node))
            logger.info("  Cores per GPU: {}".format(cores_per_gpu))
            logger.info("  Memory: {} GB/node".format(memory_gb))
            
            return config
            
        except Exception as e:
            logger.error("GPU node configuration detection failed: {}".format(e))
            return None
    
    def _parse_gres_gpu(self, partition_info):
        """Parse GPU count from partition GRES information."""
        patterns = [
            r'gpu:(\d+)',
            r'gpu:\w+:(\d+)',
            r'gres/gpu=(\d+)',
        ]
        
        for pattern in patterns:
            match = re.search(pattern, partition_info, re.IGNORECASE)
            if match:
                return int(match.group(1))
        
        return 0
    
    def _parse_cpus_per_node(self, partition_info):
        """Parse CPU count from partition information."""
        patterns = [
            r'MaxCPUsPerNode=(\d+)',
            r'CPUs=(\d+)',
        ]
        
        for pattern in patterns:
            match = re.search(pattern, partition_info)
            if match:
                return int(match.group(1))
        
        return 32  # Default
    
    def _get_sample_node(self, partition):
        """Get a sample node name from partition."""
        try:
            cmd = ["sinfo", "-p", partition, "-N", "-h", "-o", "%n"]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
            
            if result.returncode == 0 and result.stdout.strip():
                nodes = result.stdout.strip().split('\n')
                return nodes[0] if nodes else None
        except Exception:
            pass
        
        return None
    
    def _query_node_resources(self, node_name):
        """Query CPU and memory resources from a specific node."""
        try:
            cmd = ["scontrol", "show", "node", node_name]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode != 0:
                return 32, 0
            
            node_info = result.stdout
            
            # Parse CPUs
            cpus = 32
            cpu_match = re.search(r'CPUTot=(\d+)', node_info)
            if cpu_match:
                cpus = int(cpu_match.group(1))
            
            # Parse memory (in MB)
            memory_mb = 0
            mem_match = re.search(r'RealMemory=(\d+)', node_info)
            if mem_match:
                memory_mb = int(mem_match.group(1))
            
            return cpus, memory_mb // 1024
            
        except Exception:
            return 32, 0
    
    def _detect_gpu_model(self, node_name):
        """Detect GPU model and vendor."""
        # Try nvidia-smi first
        try:
            cmd = ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
            
            if result.returncode == 0 and result.stdout.strip():
                gpu_name = result.stdout.strip().split('\n')[0]
                return gpu_name, 'nvidia'
        except Exception:
            pass
        
        # Try rocm-smi
        try:
            cmd = ["rocm-smi", "--showproductname"]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
            
            if result.returncode == 0:
                for line in result.stdout.split('\n'):
                    if 'Card series' in line or 'Card model' in line:
                        parts = line.split(':')
                        if len(parts) > 1:
                            return parts[1].strip(), 'amd'
        except Exception:
            pass
        
        return "Unknown", "unknown"
    
    def generate_slurm_directives(self, num_nodes, use_all_cores=True):
        """
        Generate SLURM directives for GPU jobs.
        
        Key insight: --ntasks-per-node allocates CPUs, not MPI tasks!
        """
        if not self.node_config:
            raise ValueError("Node configuration not detected. Call detect_gpu_node_config() first.")
        
        config = self.node_config
        
        directives = {
            'nodes': str(num_nodes),
            'partition': config.partition,
            'ntasks-per-node': str(config.cpus_per_node if use_all_cores else config.gpus_per_node),
            'cpus-per-task': '1',
            'gres': 'gpu:{}'.format(config.gpus_per_node),
        }
        
        logger.info("✓ SLURM directives generated:")
        logger.info("  --nodes={}".format(num_nodes))
        logger.info("  --ntasks-per-node={} (CPU cores)".format(directives['ntasks-per-node']))
        logger.info("  --gres=gpu:{}".format(config.gpus_per_node))
        logger.info("  → Actual MPI tasks: {} per node".format(config.gpus_per_node))
        logger.info("  → Each task gets {} cores".format(config.cores_per_gpu))
        
        return directives
    
    def generate_mpirun_command(self, executable, args="", binding_script=None, verbose=False):
        """
        Generate mpirun command for GPU jobs.
        
        Automatically detects MPI implementation and uses appropriate syntax:
        - OpenMPI: Uses --map-by ppr:N:node --bind-to core --cpus-per-proc M
        - Intel MPI: Uses similar syntax WITHOUT --report-bindings
        - Other: Safe defaults
        
        Args:
            executable: Path to executable
            args: Command-line arguments
            binding_script: Optional binding script
            verbose: Add verbose/debugging flags (if supported by MPI)
        """
        if not self.node_config:
            raise ValueError("Node configuration not detected.")
        
        config = self.node_config
        
        # Import MPI detector (conditional for compatibility)
        try:
            from utils.mpi_detector import MPIDetector
            HAS_MPI_DETECTOR = True
        except ImportError:
            HAS_MPI_DETECTOR = False
        
        # Detect MPI implementation
        mpi_info = None
        if HAS_MPI_DETECTOR:
            try:
                detector = MPIDetector()
                mpi_info = detector.detect()
            except Exception as e:
                logger.warning("MPI detection failed: {}, using safe defaults".format(e))
        
        # Build command
        cmd_parts = []
        
        # Launcher name
        if mpi_info:
            cmd_parts.append(mpi_info.launcher_name)
        else:
            cmd_parts.append("mpirun")
        
        # Number of processes
        cmd_parts.append("-np {}".format(config.gpus_per_node))
        
        # Process mapping (if supported) - FIXED: use separate mapping and binding
        if mpi_info and mpi_info.supports_map_by:
            cmd_parts.append("--map-by ppr:{}:node".format(config.gpus_per_node))
            cmd_parts.append("--bind-to core")
            if config.cores_per_gpu > 1:
                cmd_parts.append("--cpus-per-proc {}".format(config.cores_per_gpu))
            logger.info("  Using GPU-aware mapping: --map-by ppr:{}:node --bind-to core".format(
                config.gpus_per_node))
        else:
            logger.info("  MPI-agnostic mode: relying on SLURM for process placement")
        
        # Binding reports (if supported and requested)
        if verbose:
            if mpi_info and mpi_info.supports_report_bindings:
                cmd_parts.append("--report-bindings")
                logger.info("  ✓ Added --report-bindings (supported by {} MPI)".format(
                    mpi_info.vendor
                ))
            else:
                mpi_vendor = mpi_info.vendor if mpi_info else 'unknown'
                logger.info("  ✗ Skipping --report-bindings (not supported by {} MPI)".format(
                    mpi_vendor
                ))
        
        # Binding script (if provided)
        if binding_script:
            cmd_parts.append(binding_script)
        
        # Executable
        cmd_parts.append(executable)
        
        # Arguments
        if args:
            cmd_parts.append(args)
        
        return ' '.join(cmd_parts)
    
    def generate_srun_command(self, executable, args=""):
        """
        Generate srun command for GPU jobs.
        """
        if not self.node_config:
            raise ValueError("Node configuration not detected.")
        
        config = self.node_config
        
        cmd_parts = [
            "srun",
            "--ntasks-per-node={}".format(config.gpus_per_node),
            "--gpus-per-node={}".format(config.gpus_per_node),
            "--gpu-bind=closest",
            "--cpu-bind=cores",
            executable,
        ]
        
        if args:
            cmd_parts.append(args)
        
        return ' '.join(cmd_parts)
    
    def generate_binding_script(self, output_path="./gpu_bind.sh"):
        """
        Generate GPU binding script similar to user's bind.sh.
        """
        if not self.node_config:
            raise ValueError("Node configuration not detected.")
        
        config = self.node_config
        
        if config.gpu_vendor == 'nvidia':
            device_var = 'CUDA_VISIBLE_DEVICES'
        elif config.gpu_vendor == 'amd':
            device_var = 'ROCR_VISIBLE_DEVICES'
        else:
            device_var = 'GPU_DEVICE_ORDINAL'
        
        script = '''#!/bin/bash
# GPU Binding Script
# Auto-generated by HPC-ScaleTest Advanced GPU Manager
#
# Maps MPI ranks to GPUs (1:1)
# System: {} with {} GPUs per node

# Get local rank (rank within node)
if [ -n "$SLURM_LOCALID" ]; then
    LOCAL_RANK=$SLURM_LOCALID
elif [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
    LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ -n "$MPI_LOCALRANKID" ]; then
    LOCAL_RANK=$MPI_LOCALRANKID
else
    echo "Warning: Could not determine local rank, using 0"
    LOCAL_RANK=0
fi

# Set GPU device based on local rank
export {}=$LOCAL_RANK

# Optional: Report binding
echo "Rank $SLURM_PROCID (local rank $LOCAL_RANK) → GPU $LOCAL_RANK on $(hostname)"

# Execute the actual application
exec "$@"
'''.format(config.partition, config.gpus_per_node, device_var)
        
        with open(output_path, 'w') as f:
            f.write(script)
        
        import os
        os.chmod(output_path, 0o755)
        
        logger.info("✓ Generated GPU binding script: {}".format(output_path))
        
        return output_path
    
    def get_optimal_decomposition(self, scaling_type='weak', dimensions=2):
        """Get optimal MPI decomposition for GPU jobs."""
        if not self.node_config:
            raise ValueError("Node configuration not detected.")
        
        gpus = self.node_config.gpus_per_node
        
        if dimensions == 2:
            if gpus == 4:
                return [2, 2, 1]
            elif gpus == 8:
                return [4, 2, 1]
            elif gpus == 6:
                return [3, 2, 1]
            else:
                return [gpus, 1, 1]
        elif dimensions == 3:
            if gpus == 8:
                return [2, 2, 2]
            elif gpus == 4:
                return [2, 2, 1]
            else:
                return [gpus, 1, 1]
        
        return [gpus, 1, 1]


def test_gpu_manager():
    """Test the GPU manager."""
    import sys
    
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    
    if len(sys.argv) < 2:
        print("Usage: python3 advanced_gpu_manager_py36.py <partition>")
        print("Example: python3 advanced_gpu_manager_py36.py booster")
        print("         python3 advanced_gpu_manager_py36.py boost_usr_prod")
        sys.exit(1)
    
    partition = sys.argv[1]
    
    print("="*70)
    print("Advanced GPU Resource Manager - Test")
    print("="*70)
    print()
    
    manager = AdvancedGPUManager()
    
    # Detect configuration
    config = manager.detect_gpu_node_config(partition)
    
    if not config:
        print("\n❌ Could not detect GPU configuration for partition '{}'".format(partition))
        sys.exit(1)
    
    print()
    print("="*70)
    print("SLURM Directives for 4 nodes:")
    print("="*70)
    directives = manager.generate_slurm_directives(num_nodes=4)
    for key, value in directives.items():
        print("  #SBATCH --{}={}".format(key, value))
    
    print()
    print("="*70)
    print("MPI Launch Commands (mpirun):")
    print("="*70)
    print("  Production mode (no verbose output):")
    mpirun_cmd = manager.generate_mpirun_command(
        executable="$BINARY/app",
        args="input.txt",
        verbose=False
    )
    print("    {}".format(mpirun_cmd))
    
    print()
    print("  Debug mode (verbose, if supported):")
    mpirun_cmd_verbose = manager.generate_mpirun_command(
        executable="$BINARY/app",
        args="input.txt",
        verbose=True
    )
    print("    {}".format(mpirun_cmd_verbose))
    
    print()
    print("="*70)
    print("MPI Launch Command (srun):")
    print("="*70)
    srun_cmd = manager.generate_srun_command(
        executable="$BINARY/app",
        args="input.txt"
    )
    print("  {}".format(srun_cmd))
    
    print()
    print("="*70)
    print("Optimal MPI Decompositions:")
    print("="*70)
    decomp_2d = manager.get_optimal_decomposition(scaling_type='weak', dimensions=2)
    decomp_3d = manager.get_optimal_decomposition(scaling_type='weak', dimensions=3)
    print("  2D weak scaling: {} = {} tasks".format(decomp_2d, decomp_2d[0]*decomp_2d[1]*decomp_2d[2]))
    print("  3D weak scaling: {} = {} tasks".format(decomp_3d, decomp_3d[0]*decomp_3d[1]*decomp_3d[2]))
    
    print()
    print("="*70)
    print("Generating GPU binding script...")
    print("="*70)
    script_path = manager.generate_binding_script("./gpu_bind_test.sh")
    print("  Created: {}".format(script_path))
    print()


if __name__ == '__main__':
    test_gpu_manager()
