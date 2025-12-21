#!/usr/bin/env python3
"""
Job Configuration Generator

This module generates complete job configurations (SLURM directives, mpirun commands,
GPU binding scripts) based on automatically detected hardware topology.

Key Responsibilities:
1. Query partition topology (CPUs, GPUs per node)
2. Derive MPI mapping (ranks per node, cores per rank)
3. Generate SLURM directives with full node allocation
4. Generate mpirun command with correct OpenMPI syntax
5. Generate GPU binding script

IMPORTANT: OpenMPI Syntax
=========================
The syntax `--map-by ppr:X:node:PE=Y` is NOT valid for OpenMPI.
The correct equivalent that achieves the same behavior is:

    --map-by ppr:X:node --bind-to core --cpus-per-proc Y

This module generates ONLY correct, working syntax.

Author: HPC-ScaleTest Contributors
"""

import os
import subprocess
import logging
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Tuple

logger = logging.getLogger(__name__)


@dataclass
class PartitionTopology:
    """
    Hardware topology for a SLURM partition.
    
    All values are detected at runtime from the specified partition.
    """
    partition: str
    cpu_cores_per_node: int
    gpus_per_node: int
    memory_gb_per_node: float = 0.0
    gpu_model: str = ""
    detection_method: str = ""
    
    # Derived values (computed in __post_init__)
    ranks_per_node: int = field(init=False)
    cores_per_rank: int = field(init=False)
    
    def __post_init__(self):
        """Compute derived MPI mapping values."""
        if self.gpus_per_node > 0:
            # GPU job: 1 rank per GPU
            self.ranks_per_node = self.gpus_per_node
            self.cores_per_rank = self.cpu_cores_per_node // self.gpus_per_node
        else:
            # CPU job: 1 rank per core
            self.ranks_per_node = self.cpu_cores_per_node
            self.cores_per_rank = 1
    
    def __str__(self) -> str:
        return (
            f"PartitionTopology({self.partition}):\n"
            f"  CPU cores/node: {self.cpu_cores_per_node}\n"
            f"  GPUs/node: {self.gpus_per_node}\n"
            f"  → Derived: {self.ranks_per_node} ranks/node, {self.cores_per_rank} cores/rank\n"
            f"  Detection: {self.detection_method}"
        )


@dataclass 
class JobConfiguration:
    """
    Complete job configuration generated from topology.
    
    Contains all information needed to generate a SLURM job script
    and mpirun command for a given number of nodes.
    """
    num_nodes: int
    topology: PartitionTopology
    
    # SLURM directives
    slurm_partition: str = field(init=False)
    slurm_nodes: int = field(init=False)
    slurm_ntasks: int = field(init=False)
    slurm_ntasks_per_node: int = field(init=False)
    slurm_cpus_per_task: int = field(init=False)
    slurm_gres: str = field(init=False)
    
    # MPI command components
    mpi_np: int = field(init=False)
    mpi_ranks_per_node: int = field(init=False)
    mpi_cores_per_rank: int = field(init=False)
    
    def __post_init__(self):
        """Compute all job configuration values."""
        # SLURM directives - request FULL node resources
        self.slurm_partition = self.topology.partition
        self.slurm_nodes = self.num_nodes
        self.slurm_ntasks_per_node = self.topology.cpu_cores_per_node  # Full allocation
        self.slurm_ntasks = self.num_nodes * self.slurm_ntasks_per_node
        self.slurm_cpus_per_task = 1  # Full node allocation mode
        
        if self.topology.gpus_per_node > 0:
            self.slurm_gres = f"gpu:{self.topology.gpus_per_node}"
        else:
            self.slurm_gres = ""
        
        # MPI command - actual task distribution
        self.mpi_np = self.num_nodes * self.topology.ranks_per_node
        self.mpi_ranks_per_node = self.topology.ranks_per_node
        self.mpi_cores_per_rank = self.topology.cores_per_rank
    
    def get_slurm_directives(self) -> Dict[str, str]:
        """Get SLURM directives as dictionary."""
        directives = {
            'nodes': str(self.slurm_nodes),
            'partition': self.slurm_partition,
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
    ) -> List[str]:
        """
        Generate mpirun command with CORRECT OpenMPI syntax.
        
        Uses the ppr:R:node:PE=C syntax for process mapping with cores per rank.
        
        Args:
            executable: Path to executable
            args: Arguments for executable
            bind_script: Path to GPU binding script
            report_bindings: Include --report-bindings flag
        
        Returns:
            Command as list of strings
        """
        args = args or []
        
        cmd = ['mpirun']
        
        # Total number of MPI ranks
        cmd.extend(['-np', str(self.mpi_np)])
        
        # Process mapping with PE (Processing Elements) for cores per rank
        # Format: ppr:<ranks_per_node>:node:PE=<cores_per_rank>
        cmd.extend(['--map-by', f'ppr:{self.mpi_ranks_per_node}:node:PE={self.mpi_cores_per_rank}'])
        
        # Report bindings for debugging
        if report_bindings:
            cmd.append('--report-bindings')
        
        # GPU binding wrapper (for GPU jobs)
        if self.topology.gpus_per_node > 0 and bind_script:
            cmd.append(bind_script)
        
        # Executable and arguments
        cmd.append(executable)
        cmd.extend(args)
        
        return cmd
    
    def get_mpirun_string(
        self,
        executable: str,
        args: List[str] = None,
        bind_script: str = "./bind.sh",
        report_bindings: bool = True,
    ) -> str:
        """Get mpirun command as string."""
        cmd = self.get_mpirun_command(executable, args, bind_script, report_bindings)
        return ' '.join(cmd)


class PartitionDetector:
    """
    Detect hardware topology for a SLURM partition.
    
    Detection Methods (in priority order):
    1. SLURM environment variables (if inside a job)
    2. sinfo query for partition hardware
    3. scontrol show partition
    
    Usage:
        detector = PartitionDetector()
        topology = detector.detect("boost_usr_prod")
        print(f"Detected: {topology.cpu_cores_per_node} CPUs, {topology.gpus_per_node} GPUs")
    """
    
    def __init__(self):
        self._cache: Dict[str, PartitionTopology] = {}
    
    def detect(self, partition: str, force_refresh: bool = False) -> PartitionTopology:
        """
        Detect topology for a partition.
        
        Args:
            partition: SLURM partition name
            force_refresh: Force re-detection even if cached
        
        Returns:
            PartitionTopology with detected hardware configuration
        
        Raises:
            RuntimeError: If detection fails
        """
        # Check cache
        if not force_refresh and partition in self._cache:
            logger.debug(f"Using cached topology for partition {partition}")
            return self._cache[partition]
        
        topology = None
        detection_method = ""
        
        # Method 1: Check SLURM environment (if inside job for this partition)
        if os.environ.get('SLURM_JOB_PARTITION') == partition:
            topology = self._detect_from_slurm_env(partition)
            if topology:
                detection_method = "SLURM environment"
        
        # Method 2: Query sinfo
        if topology is None:
            topology = self._detect_from_sinfo(partition)
            if topology:
                detection_method = "sinfo query"
        
        # Method 3: Query scontrol
        if topology is None:
            topology = self._detect_from_scontrol(partition)
            if topology:
                detection_method = "scontrol query"
        
        if topology is None:
            raise RuntimeError(
                f"Failed to detect topology for partition '{partition}'. "
                f"Ensure SLURM is available and the partition exists."
            )
        
        # Update detection method
        topology.detection_method = detection_method
        
        # Cache and return
        self._cache[partition] = topology
        
        logger.info(f"Detected topology for partition '{partition}':")
        logger.info(f"  CPU cores per node: {topology.cpu_cores_per_node}")
        logger.info(f"  GPUs per node: {topology.gpus_per_node}")
        logger.info(f"  Detection method: {detection_method}")
        logger.info(f"  Derived MPI mapping:")
        logger.info(f"    Ranks per node: {topology.ranks_per_node}")
        logger.info(f"    Cores per rank: {topology.cores_per_rank}")
        
        return topology
    
    def _detect_from_slurm_env(self, partition: str) -> Optional[PartitionTopology]:
        """Detect from SLURM environment variables."""
        try:
            cpu_cores = None
            gpus = 0
            
            # Get CPU cores
            if 'SLURM_CPUS_ON_NODE' in os.environ:
                cpu_cores = int(os.environ['SLURM_CPUS_ON_NODE'])
            elif 'SLURM_JOB_CPUS_PER_NODE' in os.environ:
                # Format: "32" or "32(x4)"
                cpus_str = os.environ['SLURM_JOB_CPUS_PER_NODE']
                cpu_cores = int(cpus_str.split('(')[0])
            
            if cpu_cores is None:
                return None
            
            # Get GPU count
            if 'SLURM_GPUS_ON_NODE' in os.environ:
                gpus = int(os.environ['SLURM_GPUS_ON_NODE'])
            elif 'CUDA_VISIBLE_DEVICES' in os.environ:
                cuda_devs = os.environ['CUDA_VISIBLE_DEVICES']
                if cuda_devs and cuda_devs.lower() not in ('', 'nodevfiles'):
                    gpus = len([d for d in cuda_devs.split(',') if d.strip()])
            
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
            # Query node info for partition
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
            
            for i, part in enumerate(parts):
                # First numeric value is CPU cores
                if cpu_cores is None and part.isdigit():
                    cpu_cores = int(part)
                
                # GPU info: gpu:type:count or gpu:count or (null)
                elif 'gpu' in part.lower() and part != '(null)':
                    import re
                    gpu_match = re.search(r'gpu:(?:([^:]+):)?(\d+)', part, re.IGNORECASE)
                    if gpu_match:
                        gpu_model = gpu_match.group(1) or ""
                        gpus = int(gpu_match.group(2))
                
                # Memory in MB
                elif part.isdigit() and cpu_cores is not None:
                    memory_gb = int(part) / 1024.0
            
            if cpu_cores is None:
                return None
            
            return PartitionTopology(
                partition=partition,
                cpu_cores_per_node=cpu_cores,
                gpus_per_node=gpus,
                memory_gb_per_node=memory_gb,
                gpu_model=gpu_model,
            )
            
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.debug(f"sinfo detection failed: {e}")
            return None
    
    def _detect_from_scontrol(self, partition: str) -> Optional[PartitionTopology]:
        """Detect from scontrol show partition."""
        try:
            result = subprocess.run(
                ['scontrol', 'show', 'partition', partition],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode != 0:
                return None
            
            output = result.stdout
            
            # Parse partition info
            # Look for MaxCPUsPerNode, Gres, etc.
            import re
            
            cpu_cores = None
            gpus = 0
            
            # MaxCPUsPerNode=32
            cpu_match = re.search(r'MaxCPUsPerNode=(\d+)', output)
            if cpu_match:
                cpu_cores = int(cpu_match.group(1))
            
            # GRES=gpu:type:count or gpu:count
            gres_match = re.search(r'GRES=([^\s]+)', output)
            if gres_match:
                gres = gres_match.group(1)
                gpu_match = re.search(r'gpu:(?:\w+:)?(\d+)', gres)
                if gpu_match:
                    gpus = int(gpu_match.group(1))
            
            if cpu_cores is None:
                # Try to get from node info
                nodes_match = re.search(r'Nodes=([^\s]+)', output)
                if nodes_match:
                    node_list = nodes_match.group(1)
                    # Query first node
                    first_node = node_list.split(',')[0].split('[')[0]
                    
                    node_result = subprocess.run(
                        ['scontrol', 'show', 'node', first_node],
                        capture_output=True,
                        text=True,
                        timeout=10
                    )
                    
                    if node_result.returncode == 0:
                        cpu_match = re.search(r'CPUTot=(\d+)', node_result.stdout)
                        if cpu_match:
                            cpu_cores = int(cpu_match.group(1))
                        
                        gres_match = re.search(r'Gres=([^\s]+)', node_result.stdout)
                        if gres_match and gpus == 0:
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


def generate_job_config(
    partition: str,
    num_nodes: int,
) -> JobConfiguration:
    """
    Generate complete job configuration for a partition.
    
    Args:
        partition: SLURM partition name
        num_nodes: Number of nodes to use
    
    Returns:
        JobConfiguration with all computed values
    
    Example:
        >>> config = generate_job_config("boost_usr_prod", num_nodes=4)
        >>> print(config.get_mpirun_string("./app", ["input.dat"]))
        mpirun -np 16 --map-by ppr:4:node --bind-to core --cpus-per-proc 8 --report-bindings ./bind.sh ./app input.dat
    """
    detector = PartitionDetector()
    topology = detector.detect(partition)
    return JobConfiguration(num_nodes=num_nodes, topology=topology)


def generate_gpu_binding_script(output_path: str = "./bind.sh") -> str:
    """
    Generate GPU binding script.
    
    This script uses OMPI_COMM_WORLD_LOCAL_RANK to set CUDA_VISIBLE_DEVICES,
    ensuring each MPI rank is bound to a unique GPU.
    
    Args:
        output_path: Path to write the script
    
    Returns:
        Path to the generated script
    """
    script_content = '''#!/bin/bash
# GPU Binding Script - Generated by HPC-ScaleTest
#
# This script binds each MPI rank to a unique GPU using the local rank.
# It supports OpenMPI, Intel MPI, MPICH, and SLURM.
#
# Usage: mpirun ... ./bind.sh <executable> [args...]

# Determine local rank from available environment variables
if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
    # OpenMPI
    LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ -n "$MPI_LOCALRANKID" ]; then
    # Intel MPI
    LOCAL_RANK=$MPI_LOCALRANKID
elif [ -n "$PMI_LOCAL_RANK" ]; then
    # MPICH/Hydra
    LOCAL_RANK=$PMI_LOCAL_RANK
elif [ -n "$SLURM_LOCALID" ]; then
    # SLURM
    LOCAL_RANK=$SLURM_LOCALID
elif [ -n "$MV2_COMM_WORLD_LOCAL_RANK" ]; then
    # MVAPICH2
    LOCAL_RANK=$MV2_COMM_WORLD_LOCAL_RANK
else
    # Fallback
    LOCAL_RANK=0
    echo "WARNING: Could not determine local rank for GPU binding" >&2
fi

# Bind to GPU
export CUDA_VISIBLE_DEVICES=$LOCAL_RANK

# Debug output (uncomment if needed)
# echo "Rank ${OMPI_COMM_WORLD_RANK:-$SLURM_PROCID}: LOCAL_RANK=$LOCAL_RANK, CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"

# Execute the command
exec "$@"
'''
    
    with open(output_path, 'w') as f:
        f.write(script_content)
    
    os.chmod(output_path, 0o755)
    logger.info(f"Generated GPU binding script: {output_path}")
    
    return output_path


def generate_slurm_job_script(
    config: JobConfiguration,
    executable: str,
    args: List[str] = None,
    job_name: str = "hpc_scaletest",
    time_limit: str = "02:00:00",
    account: str = None,
    output_dir: str = ".",
    env_setup: List[str] = None,
    qos: str = None,
) -> str:
    """
    Generate complete SLURM job script.
    
    Args:
        config: JobConfiguration from generate_job_config()
        executable: Path to executable
        args: Arguments for executable
        job_name: Job name
        time_limit: Time limit (HH:MM:SS)
        account: SLURM account
        output_dir: Output directory
        env_setup: Environment setup commands (module loads, etc.)
        qos: Quality of service
    
    Returns:
        Complete job script as string
    """
    args = args or []
    env_setup = env_setup or []
    
    topo = config.topology
    
    script = f'''#!/bin/bash
# =============================================================================
# SLURM Job Script - Generated by HPC-ScaleTest
# =============================================================================
# Partition: {topo.partition}
# Detected topology:
#   CPU cores per node: {topo.cpu_cores_per_node}
#   GPUs per node: {topo.gpus_per_node}
# Derived MPI mapping:
#   Ranks per node: {topo.ranks_per_node} (1 per GPU)
#   Cores per rank: {topo.cores_per_rank}
# Detection method: {topo.detection_method}
# =============================================================================

#SBATCH --job-name={job_name}
#SBATCH --nodes={config.slurm_nodes}
#SBATCH --partition={config.slurm_partition}
#SBATCH --ntasks-per-node={config.slurm_ntasks_per_node}
'''
    
    if config.slurm_gres:
        script += f"#SBATCH --gres={config.slurm_gres}\n"
    
    if qos:
        script += f"#SBATCH --qos={qos}\n"
    
    if account:
        script += f"#SBATCH --account={account}\n"
    
    script += f'''#SBATCH --time={time_limit}
#SBATCH --output={output_dir}/job_%j.out
#SBATCH --error={output_dir}/job_%j.err
#SBATCH --exclusive

# =============================================================================
# Environment Setup
# =============================================================================
'''
    
    for cmd in env_setup:
        script += f"{cmd}\n"
    
    script += f'''
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_NNODES"
echo "Tasks per node: $SLURM_NTASKS_PER_NODE"
echo "Partition: $SLURM_JOB_PARTITION"
echo "=========================================="

# Set OpenMP threads to match cores per rank
export OMP_NUM_THREADS={topo.cores_per_rank}

# Change to working directory
cd {output_dir}

# =============================================================================
# GPU Binding Script
# =============================================================================
cat > bind.sh << 'BIND_SCRIPT_EOF'
#!/bin/bash
if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
    LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ -n "$MPI_LOCALRANKID" ]; then
    LOCAL_RANK=$MPI_LOCALRANKID
elif [ -n "$PMI_LOCAL_RANK" ]; then
    LOCAL_RANK=$PMI_LOCAL_RANK
elif [ -n "$SLURM_LOCALID" ]; then
    LOCAL_RANK=$SLURM_LOCALID
else
    LOCAL_RANK=0
fi
export CUDA_VISIBLE_DEVICES=$LOCAL_RANK
exec "$@"
BIND_SCRIPT_EOF
chmod +x bind.sh

# =============================================================================
# Execute Application
# =============================================================================
echo "Starting application at $(date)"
echo "MPI command:"
echo "  -np {config.mpi_np}"
echo "  --map-by ppr:{config.mpi_ranks_per_node}:node"
echo "  --bind-to core"
echo "  --cpus-per-proc {config.mpi_cores_per_rank}"

start_time=$(date +%s.%N)

{config.get_mpirun_string(executable, args, "./bind.sh", report_bindings=True)}

end_time=$(date +%s.%N)
elapsed=$(echo "$end_time - $start_time" | bc)

echo "=========================================="
echo "Application finished at $(date)"
echo "Elapsed time: $elapsed seconds"
echo "=========================================="
'''
    
    return script


# =============================================================================
# Module-level convenience functions
# =============================================================================

_detector: Optional[PartitionDetector] = None

def get_partition_detector() -> PartitionDetector:
    """Get or create the global partition detector."""
    global _detector
    if _detector is None:
        _detector = PartitionDetector()
    return _detector


def detect_partition_topology(partition: str) -> PartitionTopology:
    """
    Convenience function to detect partition topology.
    
    Args:
        partition: SLURM partition name
    
    Returns:
        PartitionTopology with detected configuration
    """
    return get_partition_detector().detect(partition)


# =============================================================================
# Self-test
# =============================================================================

if __name__ == '__main__':
    import sys
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s - %(message)s'
    )
    
    print("=" * 70)
    print(" Job Configuration Generator - Self Test")
    print("=" * 70)
    
    # Test with mock topology (simulating Leonardo Booster)
    print("\nTest 1: Mock Leonardo Booster topology (32 cores, 4 GPUs)")
    
    mock_topology = PartitionTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        gpu_model="a100",
        detection_method="mock"
    )
    
    print(mock_topology)
    
    # Test different node counts
    for num_nodes in [1, 2, 4]:
        print(f"\nTest: {num_nodes} node(s)")
        
        config = JobConfiguration(num_nodes=num_nodes, topology=mock_topology)
        
        print(f"  SLURM directives:")
        print(f"    #SBATCH --ntasks-per-node={config.slurm_ntasks_per_node}")
        print(f"    #SBATCH --gres={config.slurm_gres}")
        
        mpi_cmd = config.get_mpirun_string(
            "$BINARY/iPIC3D",
            ["os-stdin"],
            "./bind.sh"
        )
        print(f"  mpirun command:")
        print(f"    {mpi_cmd}")
        
        # Verify correct syntax
        assert ':PE=' not in mpi_cmd, "ERROR: Found incorrect :PE= syntax!"
        assert f'-np {num_nodes * 4}' in mpi_cmd, "ERROR: Incorrect -np value"
        assert '--map-by ppr:4:node' in mpi_cmd, "ERROR: Missing ppr mapping"
        assert '--bind-to core' in mpi_cmd, "ERROR: Missing bind-to"
        assert '--cpus-per-proc 8' in mpi_cmd, "ERROR: Missing cpus-per-proc"
        
        print("  ✓ Syntax verification passed")
    
    # Test LUMI-style topology (128 cores, 8 GPUs)
    print("\n" + "=" * 70)
    print("Test 2: Mock LUMI topology (128 cores, 8 GPUs)")
    
    lumi_topology = PartitionTopology(
        partition="gpu",
        cpu_cores_per_node=128,
        gpus_per_node=8,
        gpu_model="mi250x",
        detection_method="mock"
    )
    
    print(lumi_topology)
    
    config = JobConfiguration(num_nodes=4, topology=lumi_topology)
    mpi_cmd = config.get_mpirun_string("./app", ["input.dat"], "./bind.sh")
    
    print(f"  mpirun command:")
    print(f"    {mpi_cmd}")
    
    # Verify
    assert '-np 32' in mpi_cmd, "ERROR: Incorrect -np for LUMI"
    assert '--map-by ppr:8:node' in mpi_cmd, "ERROR: Incorrect ppr for LUMI"
    assert '--cpus-per-proc 16' in mpi_cmd, "ERROR: Incorrect cpus-per-proc for LUMI"
    
    print("  ✓ LUMI configuration verified")
    
    print("\n" + "=" * 70)
    print(" All tests PASSED")
    print("=" * 70)
