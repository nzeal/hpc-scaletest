"""
CPU Job Generator Module

Generates job scripts for CPU-only HPC jobs.
This module handles pure MPI jobs without GPU acceleration.

Design Principles:
==================
- NO HARDCODED VALUES: All values from configuration
- SYSTEM AGNOSTIC: Works on any HPC system with SLURM
- SIMPLE MPI: Standard mpirun/srun without GPU binding

MPI Command Format:
===================
For mpirun: mpirun -np <total_ranks> $BINARY/<executable> <input>
For srun:   srun --ntasks=<total_ranks> $BINARY/<executable> <input>

Author: HPC-ScaleTest Contributors
"""

import logging
from typing import Optional

from .base import JobGeneratorBase, JobGeneratorConfig

logger = logging.getLogger(__name__)


class CPUJobGenerator(JobGeneratorBase):
    """
    Job script generator for CPU-only jobs.
    
    For CPU jobs:
    - 1 MPI rank per CPU core (default)
    - Simple mpirun/srun without complex mapping
    - No GPU binding required
    
    Example (32 cores/node, 4 nodes):
        mpirun -np 128 $BINARY/app input.dat
    """
    
    def __init__(self, config: JobGeneratorConfig):
        """
        Initialize CPU job generator.
        
        Args:
            config: Job generator configuration
        """
        super().__init__(config)
        
        # Validate: gpus_per_node should be 0 for CPU jobs
        if config.gpus_per_node > 0:
            logger.warning(
                f"CPUJobGenerator initialized with gpus_per_node={config.gpus_per_node}. "
                f"For GPU jobs, use GPUJobGenerator instead."
            )
        
        # Set procs_per_node if not explicitly configured
        if config.procs_per_node == 0:
            config.procs_per_node = config.cpus_per_node
    
    def generate_mpi_command(
        self,
        num_nodes: int,
        total_ranks: int,
        ranks_per_node: int,
        cores_per_rank: int
    ) -> str:
        """
        Generate MPI launch command for CPU job.
        
        Args:
            num_nodes: Number of nodes
            total_ranks: Total MPI ranks
            ranks_per_node: MPI ranks per node
            cores_per_rank: CPU cores per rank (typically 1 for CPU jobs)
        
        Returns:
            MPI launch command
        """
        exe_path = f"$BINARY/{self.config.executable}" if self.config.binary_dir else self.config.executable
        input_arg = self.config.input_file if self.config.input_file else ""
        
        if self.config.launcher == "srun":
            # SLURM native launcher
            cmd = f"srun --ntasks={total_ranks} --ntasks-per-node={ranks_per_node} {exe_path} {input_arg}"
        else:
            # mpirun (OpenMPI, MPICH, etc.)
            cmd = f"mpirun -np {total_ranks} {exe_path} {input_arg}"
        
        return cmd.strip()
    
    def generate_script(
        self,
        num_nodes: int,
        job_name: str,
        working_dir: str,
        procs_decomposition: Optional[tuple] = None,
        **kwargs
    ) -> str:
        """
        Generate complete CPU job script.
        
        Args:
            num_nodes: Number of nodes
            job_name: SLURM job name
            working_dir: Working directory
            procs_decomposition: MPI decomposition (px, py, pz) if needed
            **kwargs: Additional parameters
        
        Returns:
            Complete job script as string
        """
        # Calculate MPI configuration
        ranks_per_node = self.config.procs_per_node
        total_ranks = num_nodes * ranks_per_node
        cores_per_rank = 1  # Standard for CPU-only jobs
        
        # Override total_ranks if procs_decomposition provided
        if procs_decomposition:
            total_ranks = procs_decomposition[0] * procs_decomposition[1] * procs_decomposition[2]
            ranks_per_node = total_ranks // num_nodes
        
        # Generate script components
        header = self.generate_header_comment(
            num_nodes=num_nodes,
            ranks_per_node=ranks_per_node,
            cores_per_rank=cores_per_rank,
            job_type="CPU"
        )
        
        directives = self.generate_slurm_directives(
            num_nodes=num_nodes,
            job_name=job_name,
            working_dir=working_dir,
            ntasks_per_node=self.config.cpus_per_node,  # Full allocation
            gpus_per_node=0
        )
        
        module_loads = self.generate_module_loads()
        env_exports = self.generate_env_exports()
        
        mpi_cmd = self.generate_mpi_command(
            num_nodes=num_nodes,
            total_ranks=total_ranks,
            ranks_per_node=ranks_per_node,
            cores_per_rank=cores_per_rank
        )
        
        # Build complete script
        script = header
        script += "\n".join(directives) + "\n\n"
        script += module_loads + "\n"
        script += env_exports
        
        script += f'''
# Change to working directory
cd {working_dir}

# Set binary location
export BINARY={self.config.binary_dir if self.config.binary_dir else working_dir}

echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_NNODES"
echo "Tasks per node: {ranks_per_node}"
echo "Total MPI ranks: {total_ranks}"
echo "Working directory: $PWD"
echo "BINARY: $BINARY"
echo "=========================================="

echo ""
echo "Starting: $(date)"
echo "Command: {mpi_cmd}"
echo ""

# Execute with timing
start_time="$(date -u +%s.%N)"

{mpi_cmd}

exit_code=$?
end_time="$(date -u +%s.%N)"
elapsed="$(echo "$end_time - $start_time" | bc)"

echo ""
echo "=========================================="
echo "Completed: $(date)"
echo "Exit code: $exit_code"
echo "Elapsed: $elapsed seconds"
echo "=========================================="

exit $exit_code
'''
        
        return script
