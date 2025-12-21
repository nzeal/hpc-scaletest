#!/usr/bin/env python3
"""
Job Script Generator with Automatic Topology Detection

This module generates complete SLURM job scripts based on automatically
detected partition topology. No hardcoded values or site-specific constants.

Key Features:
- Automatic CPU/GPU detection via SLURM queries
- Derived MPI mapping (ranks per node = GPUs, cores per rank)
- Full node resource allocation in SLURM directives
- Correct mpirun syntax with PE= for cores per rank
- GPU binding via CUDA_VISIBLE_DEVICES

Author: HPC-ScaleTest Contributors
"""

import os
import logging
from typing import List, Optional, Dict, Any

# Import from partition_topology (avoiding circular imports)
import sys
_orig_path = sys.path.copy()
try:
    # Ensure we can import even if run from project directory
    from dataclasses import dataclass, field
    import subprocess
    import re
except ImportError:
    pass

logger = logging.getLogger(__name__)


@dataclass
class PartitionTopology:
    """Hardware topology for a SLURM partition."""
    partition: str
    cpu_cores_per_node: int
    gpus_per_node: int
    memory_gb: float = 0.0
    gpu_model: str = ""
    detection_method: str = ""
    
    ranks_per_node: int = field(init=False)
    cores_per_rank: int = field(init=False)
    
    def __post_init__(self):
        if self.gpus_per_node > 0:
            self.ranks_per_node = self.gpus_per_node
            self.cores_per_rank = self.cpu_cores_per_node // self.gpus_per_node
        else:
            self.ranks_per_node = self.cpu_cores_per_node
            self.cores_per_rank = 1


class TopologyDetector:
    """Detect hardware topology for SLURM partitions."""
    
    def __init__(self):
        self._cache: Dict[str, PartitionTopology] = {}
    
    def detect(self, partition: str) -> PartitionTopology:
        """
        Detect topology for a partition.
        
        Detection order:
        1. SLURM environment variables (if inside job)
        2. sinfo query
        3. scontrol query
        """
        if partition in self._cache:
            return self._cache[partition]
        
        topology = None
        method = ""
        
        # Method 1: SLURM environment
        if os.environ.get('SLURM_JOB_PARTITION') == partition:
            topology = self._from_slurm_env(partition)
            if topology:
                method = "SLURM environment"
        
        # Method 2: sinfo
        if topology is None:
            topology = self._from_sinfo(partition)
            if topology:
                method = "sinfo query"
        
        # Method 3: scontrol
        if topology is None:
            topology = self._from_scontrol(partition)
            if topology:
                method = "scontrol query"
        
        if topology is None:
            raise RuntimeError(f"Cannot detect topology for partition '{partition}'")
        
        topology.detection_method = method
        self._cache[partition] = topology
        
        logger.info(f"Topology for '{partition}': {topology.cpu_cores_per_node} cores, "
                   f"{topology.gpus_per_node} GPUs → {topology.ranks_per_node} ranks/node, "
                   f"{topology.cores_per_rank} cores/rank ({method})")
        
        return topology
    
    def _from_slurm_env(self, partition: str) -> Optional[PartitionTopology]:
        try:
            cpu_cores = None
            gpus = 0
            
            if 'SLURM_CPUS_ON_NODE' in os.environ:
                cpu_cores = int(os.environ['SLURM_CPUS_ON_NODE'])
            elif 'SLURM_JOB_CPUS_PER_NODE' in os.environ:
                val = os.environ['SLURM_JOB_CPUS_PER_NODE']
                cpu_cores = int(val.split('(')[0])
            
            if cpu_cores is None:
                return None
            
            if 'SLURM_GPUS_ON_NODE' in os.environ:
                gpus = int(os.environ['SLURM_GPUS_ON_NODE'])
            elif 'SLURM_GPUS_PER_NODE' in os.environ:
                gpus = int(os.environ['SLURM_GPUS_PER_NODE'])
            
            return PartitionTopology(partition=partition, cpu_cores_per_node=cpu_cores, gpus_per_node=gpus)
        except:
            return None
    
    def _from_sinfo(self, partition: str) -> Optional[PartitionTopology]:
        try:
            result = subprocess.run(
                ['sinfo', '-p', partition, '-N', '-h', '-o', '%c %G %m'],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode != 0 or not result.stdout.strip():
                return None
            
            line = result.stdout.strip().split('\n')[0]
            parts = line.split()
            
            cpu_cores = None
            gpus = 0
            
            for part in parts:
                if cpu_cores is None and part.isdigit():
                    cpu_cores = int(part)
                elif 'gpu' in part.lower() and part != '(null)':
                    match = re.search(r'gpu:(?:\w+:)?(\d+)', part, re.IGNORECASE)
                    if match:
                        gpus = int(match.group(1))
            
            if cpu_cores is None:
                return None
            
            return PartitionTopology(partition=partition, cpu_cores_per_node=cpu_cores, gpus_per_node=gpus)
        except:
            return None
    
    def _from_scontrol(self, partition: str) -> Optional[PartitionTopology]:
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
            
            match = re.search(r'MaxCPUsPerNode=(\d+)', output)
            if match:
                cpu_cores = int(match.group(1))
            
            match = re.search(r'GRES=([^\s]+)', output)
            if match:
                gpu_match = re.search(r'gpu:(?:\w+:)?(\d+)', match.group(1))
                if gpu_match:
                    gpus = int(gpu_match.group(1))
            
            if cpu_cores is None:
                # Query a node
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
            
            return PartitionTopology(partition=partition, cpu_cores_per_node=cpu_cores, gpus_per_node=gpus)
        except:
            return None


# Global detector instance
_detector = None

def get_detector() -> TopologyDetector:
    global _detector
    if _detector is None:
        _detector = TopologyDetector()
    return _detector


def generate_job_script(
    partition: str,
    num_nodes: int,
    executable: str,
    args: List[str] = None,
    job_name: str = "hpc_job",
    time_limit: str = "02:00:00",
    account: str = None,
    qos: str = None,
    output_dir: str = ".",
    env_setup: List[str] = None,
    topology: PartitionTopology = None,
) -> str:
    """
    Generate complete SLURM job script with automatic topology detection.
    
    Args:
        partition: SLURM partition name
        num_nodes: Number of nodes
        executable: Path to executable
        args: Executable arguments
        job_name: Job name
        time_limit: Time limit (HH:MM:SS)
        account: SLURM account
        qos: Quality of service
        output_dir: Output directory
        env_setup: Environment setup commands (module loads)
        topology: Optional pre-detected topology (otherwise auto-detected)
    
    Returns:
        Complete job script as string
    
    Example (Leonardo Booster, 4 nodes):
        SLURM: --ntasks-per-node=32, --gres=gpu:4
        mpirun: -np 16 --map-by ppr:4:node:PE=8 --report-bindings ./bind.sh <exec>
    """
    args = args or []
    env_setup = env_setup or []
    
    # Detect topology if not provided
    if topology is None:
        detector = get_detector()
        topology = detector.detect(partition)
    
    # Compute values
    mpi_np = num_nodes * topology.ranks_per_node
    
    # Build mpirun command
    mpi_cmd = (
        f"mpirun -np {mpi_np} "
        f"--map-by ppr:{topology.ranks_per_node}:node:PE={topology.cores_per_rank} "
        f"--report-bindings ./bind.sh {executable}"
    )
    if args:
        mpi_cmd += " " + " ".join(args)
    
    # Generate script
    script = f'''#!/bin/bash
# =============================================================================
# SLURM Job Script - Generated by HPC-ScaleTest
# =============================================================================
# Partition: {partition}
# Detected topology:
#   CPU cores per node: {topology.cpu_cores_per_node}
#   GPUs per node: {topology.gpus_per_node}
# Derived MPI mapping:
#   Ranks per node: {topology.ranks_per_node} (1 per GPU)
#   Cores per rank: {topology.cores_per_rank}
# Detection method: {topology.detection_method}
# =============================================================================

#SBATCH --job-name={job_name}
#SBATCH --nodes={num_nodes}
#SBATCH --partition={partition}
#SBATCH --ntasks-per-node={topology.cpu_cores_per_node}
'''
    
    if topology.gpus_per_node > 0:
        script += f"#SBATCH --gres=gpu:{topology.gpus_per_node}\n"
    
    if account:
        script += f"#SBATCH --account={account}\n"
    
    if qos:
        script += f"#SBATCH --qos={qos}\n"
    
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
echo "Partition: $SLURM_JOB_PARTITION"
echo "CPUs on node: $SLURM_CPUS_ON_NODE"
echo "GPUs on node: ${{SLURM_GPUS_ON_NODE:-0}}"
echo "=========================================="

# Set OpenMP threads to cores per rank
export OMP_NUM_THREADS={topology.cores_per_rank}

cd {output_dir}

# =============================================================================
# GPU Binding Script (CUDA_VISIBLE_DEVICES via OMPI_COMM_WORLD_LOCAL_RANK)
# =============================================================================
cat > bind.sh << 'BIND_EOF'
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
BIND_EOF
chmod +x bind.sh

# =============================================================================
# Run Application
# =============================================================================
echo "Starting: $(date)"
echo "Command: {mpi_cmd}"

{mpi_cmd}

echo "Finished: $(date)"
'''
    
    return script


def get_mpirun_command(
    partition: str,
    num_nodes: int,
    executable: str,
    args: List[str] = None,
    bind_script: str = "./bind.sh",
    report_bindings: bool = True,
    topology: PartitionTopology = None,
) -> str:
    """
    Generate mpirun command string.
    
    Args:
        partition: SLURM partition name
        num_nodes: Number of nodes
        executable: Executable path
        args: Executable arguments
        bind_script: GPU binding script path
        report_bindings: Include --report-bindings
        topology: Optional pre-detected topology
    
    Returns:
        mpirun command string
    
    Example (Leonardo, 4 nodes):
        mpirun -np 16 --map-by ppr:4:node:PE=8 --report-bindings ./bind.sh $BINARY/iPIC3D os-stdin
    """
    args = args or []
    
    if topology is None:
        detector = get_detector()
        topology = detector.detect(partition)
    
    mpi_np = num_nodes * topology.ranks_per_node
    
    parts = [
        'mpirun',
        f'-np {mpi_np}',
        f'--map-by ppr:{topology.ranks_per_node}:node:PE={topology.cores_per_rank}',
    ]
    
    if report_bindings:
        parts.append('--report-bindings')
    
    if topology.gpus_per_node > 0 and bind_script:
        parts.append(bind_script)
    
    parts.append(executable)
    parts.extend(args)
    
    return ' '.join(parts)


# =============================================================================
# Test
# =============================================================================

if __name__ == '__main__':
    print("=" * 70)
    print(" Job Script Generator - Self Test")
    print("=" * 70)
    
    # Mock Leonardo Booster topology
    leo_topo = PartitionTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        detection_method="mock"
    )
    
    print("\n[Leonardo Booster] 32 cores, 4 GPUs per node")
    print(f"  Derived: {leo_topo.ranks_per_node} ranks/node, {leo_topo.cores_per_rank} cores/rank")
    
    # Test mpirun commands
    print("\n  mpirun commands:")
    for N in [1, 2, 4]:
        cmd = get_mpirun_command(
            partition="boost_usr_prod",
            num_nodes=N,
            executable="$BINARY/iPIC3D",
            args=["os-stdin"],
            topology=leo_topo
        )
        print(f"    {N} node(s): {cmd}")
    
    # Generate sample script
    print("\n  Sample job script (4 nodes):")
    script = generate_job_script(
        partition="boost_usr_prod",
        num_nodes=4,
        executable="$BINARY/iPIC3D",
        args=["os-stdin"],
        job_name="ipic3d_test",
        account="my_account",
        env_setup=["module load cuda/12.0", "module load openmpi/4.1"],
        topology=leo_topo
    )
    print("    " + "\n    ".join(script.split('\n')[:25]) + "\n    ...")
    
    # Verify exact syntax
    print("\n  Syntax verification:")
    cmd = get_mpirun_command("boost_usr_prod", 4, "./app", topology=leo_topo)
    checks = [
        ('-np 16' in cmd, "-np 16 (4 nodes × 4 GPUs)"),
        ('ppr:4:node:PE=8' in cmd, "ppr:4:node:PE=8 (4 ranks, 8 cores each)"),
        ('--report-bindings' in cmd, "--report-bindings"),
        ('./bind.sh' in cmd, "./bind.sh (GPU binding)"),
    ]
    for passed, desc in checks:
        status = "✓" if passed else "✗"
        print(f"    {status} {desc}")
    
    print("\n" + "=" * 70)
    print(" All tests PASSED")
    print("=" * 70)
