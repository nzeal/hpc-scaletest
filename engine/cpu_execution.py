#!/usr/bin/env python3
"""
engine/cpu_execution.py - CPU-Only Execution Module

Handles execution logic for CPU-only HPC jobs (no GPUs).

For CPU-only jobs:
- One MPI rank per CPU core (typical)
- No GPU binding required
- Standard mpirun command generation

Author: HPC-ScaleTest Contributors
"""

import logging
from dataclasses import dataclass
from typing import List, Optional, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class CPUJobConfig:
    """Configuration for a CPU-only job."""
    
    # Node configuration
    num_nodes: int
    cpus_per_node: int
    
    # MPI configuration (derived)
    mpi_ranks_per_node: int = 0
    cores_per_rank: int = 1
    total_mpi_ranks: int = 0
    
    # Problem decomposition
    procs_x: int = 1
    procs_y: int = 1
    procs_z: int = 1
    
    def __post_init__(self):
        """Compute derived values."""
        if self.mpi_ranks_per_node == 0:
            # Default: one rank per core
            self.mpi_ranks_per_node = self.cpus_per_node
        
        if self.total_mpi_ranks == 0:
            self.total_mpi_ranks = self.num_nodes * self.mpi_ranks_per_node
        
        # Verify consistency
        expected_total = self.procs_x * self.procs_y * self.procs_z
        if expected_total != self.total_mpi_ranks:
            logger.warning(
                f"MPI decomposition {self.procs_x}×{self.procs_y}×{self.procs_z}={expected_total} "
                f"doesn't match total_mpi_ranks={self.total_mpi_ranks}"
            )


class CPUExecutionEngine:
    """
    Execution engine for CPU-only jobs.
    
    This handles:
    - Job script generation for CPU partitions
    - MPI command generation without GPU binding
    - Resource allocation for CPU-only nodes
    """
    
    def __init__(self, partition: str = ""):
        """
        Initialize CPU execution engine.
        
        Args:
            partition: SLURM partition name
        """
        self.partition = partition
        self._cpus_per_node: Optional[int] = None
    
    def detect_cpus_per_node(self) -> int:
        """
        Detect CPUs per node from SLURM or system.
        
        Returns:
            Number of CPU cores per node
        """
        import os
        import subprocess
        
        # Method 1: SLURM environment
        if 'SLURM_CPUS_ON_NODE' in os.environ:
            self._cpus_per_node = int(os.environ['SLURM_CPUS_ON_NODE'])
            logger.info(f"Detected CPUs from SLURM: {self._cpus_per_node}")
            return self._cpus_per_node
        
        # Method 2: sinfo query
        if self.partition:
            try:
                result = subprocess.run(
                    ['sinfo', '-p', self.partition, '-N', '-h', '-o', '%c'],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                if result.returncode == 0 and result.stdout.strip():
                    self._cpus_per_node = int(result.stdout.strip().split()[0])
                    logger.info(f"Detected CPUs from sinfo: {self._cpus_per_node}")
                    return self._cpus_per_node
            except Exception as e:
                logger.debug(f"sinfo detection failed: {e}")
        
        # Method 3: System introspection
        self._cpus_per_node = os.cpu_count() or 1
        logger.info(f"Detected CPUs from system: {self._cpus_per_node}")
        return self._cpus_per_node
    
    @property
    def cpus_per_node(self) -> int:
        """Get cached or detect CPUs per node."""
        if self._cpus_per_node is None:
            self.detect_cpus_per_node()
        return self._cpus_per_node
    
    def create_job_config(
        self,
        num_nodes: int,
        procs_decomp: Tuple[int, int, int] = None,
        ranks_per_node: int = None,
        cores_per_rank: int = 1
    ) -> CPUJobConfig:
        """
        Create job configuration for CPU-only execution.
        
        Args:
            num_nodes: Number of nodes
            procs_decomp: MPI decomposition (px, py, pz)
            ranks_per_node: Override ranks per node
            cores_per_rank: Cores per MPI rank (for hybrid MPI+OpenMP)
            
        Returns:
            CPUJobConfig with all values computed
        """
        cpus = self.cpus_per_node
        
        # Default: one rank per core
        if ranks_per_node is None:
            ranks_per_node = cpus // cores_per_rank
        
        total_ranks = num_nodes * ranks_per_node
        
        # Compute or validate decomposition
        if procs_decomp is None:
            procs_decomp = self._compute_decomposition(total_ranks)
        
        return CPUJobConfig(
            num_nodes=num_nodes,
            cpus_per_node=cpus,
            mpi_ranks_per_node=ranks_per_node,
            cores_per_rank=cores_per_rank,
            total_mpi_ranks=total_ranks,
            procs_x=procs_decomp[0],
            procs_y=procs_decomp[1],
            procs_z=procs_decomp[2]
        )
    
    def _compute_decomposition(self, total_ranks: int) -> Tuple[int, int, int]:
        """
        Compute a balanced MPI decomposition.
        
        Args:
            total_ranks: Total MPI processes
            
        Returns:
            Tuple (px, py, pz) decomposition
        """
        import math
        
        # Try to find balanced factorization
        factors = []
        n = total_ranks
        for p in [2, 3, 5, 7]:
            while n % p == 0:
                factors.append(p)
                n //= p
        if n > 1:
            factors.append(n)
        
        # Distribute factors into 3 dimensions
        px, py, pz = 1, 1, 1
        for i, f in enumerate(sorted(factors, reverse=True)):
            if i % 3 == 0:
                px *= f
            elif i % 3 == 1:
                py *= f
            else:
                pz *= f
        
        return (px, py, pz)
    
    def generate_mpirun_command(
        self,
        config: CPUJobConfig,
        executable: str,
        args: List[str] = None,
        report_bindings: bool = True
    ) -> List[str]:
        """
        Generate mpirun command for CPU-only job.
        
        Args:
            config: CPU job configuration
            executable: Path to executable
            args: Arguments for executable
            report_bindings: Include --report-bindings
            
        Returns:
            List of command components
        """
        args = args or []
        
        cmd = [
            'mpirun',
            '-np', str(config.total_mpi_ranks),
            '--map-by', f'ppr:{config.mpi_ranks_per_node}:node:PE={config.cores_per_rank}',
        ]
        
        if report_bindings:
            cmd.append('--report-bindings')
        
        cmd.append(executable)
        cmd.extend(args)
        
        return cmd
    
    def generate_mpirun_string(
        self,
        config: CPUJobConfig,
        executable: str,
        args: List[str] = None,
        report_bindings: bool = True
    ) -> str:
        """Generate mpirun command as string."""
        cmd = self.generate_mpirun_command(config, executable, args, report_bindings)
        return ' '.join(cmd)
    
    def generate_slurm_directives(
        self,
        config: CPUJobConfig,
        job_name: str,
        time_limit: str = "02:00:00",
        account: str = None,
        qos: str = None,
        exclusive: bool = True
    ) -> List[str]:
        """
        Generate SLURM directives for CPU-only job.
        
        Args:
            config: CPU job configuration
            job_name: SLURM job name
            time_limit: Time limit (HH:MM:SS)
            account: SLURM account
            qos: Quality of service
            exclusive: Exclusive node access
            
        Returns:
            List of SLURM directive lines
        """
        directives = [
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --nodes={config.num_nodes}",
            f"#SBATCH --ntasks-per-node={config.cpus_per_node}",  # Full node
            f"#SBATCH --time={time_limit}",
        ]
        
        if self.partition:
            directives.append(f"#SBATCH --partition={self.partition}")
        
        if account:
            directives.append(f"#SBATCH --account={account}")
        
        if qos:
            directives.append(f"#SBATCH --qos={qos}")
        
        if exclusive:
            directives.append("#SBATCH --exclusive")
        
        return directives


# =============================================================================
# Self-test
# =============================================================================

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
    
    print("=" * 70)
    print(" CPU Execution Module - Self Test")
    print("=" * 70)
    
    engine = CPUExecutionEngine(partition="cpu_partition")
    engine._cpus_per_node = 64  # Mock detection
    
    print(f"\n  Partition: cpu_partition")
    print(f"  CPUs per node: {engine.cpus_per_node}")
    
    # Test configurations
    for nodes in [1, 2, 4]:
        config = engine.create_job_config(num_nodes=nodes)
        cmd = engine.generate_mpirun_string(config, "./app", ["input.dat"])
        
        print(f"\n  {nodes} node(s):")
        print(f"    Total ranks: {config.total_mpi_ranks}")
        print(f"    Ranks/node: {config.mpi_ranks_per_node}")
        print(f"    Command: {cmd}")
        
        assert config.total_mpi_ranks == nodes * 64
        assert f"-np {nodes * 64}" in cmd
    
    print("\n" + "=" * 70)
    print(" All tests PASSED")
    print("=" * 70)
