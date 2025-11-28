"""
Job script builder for HPC schedulers.

Extracts job script generation logic from runner.py into a reusable module.
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional, List
from core.config import JobConfig  # Fixed import

logger = logging.getLogger(__name__)

# Try to import SlurmDetector for partition auto-selection
try:
    from utils.slurm_detector import get_slurm_detector
    HAS_SLURM_DETECTOR = True
except ImportError:
    HAS_SLURM_DETECTOR = False


class JobScriptBuilder:
    """Builds job submission scripts for HPC schedulers."""
    
    def __init__(
        self,
        scheduler: str = "slurm",
        partition: Optional[str] = None,
        account: Optional[str] = None,
        time_limit: str = "01:00:00",
        modules: Optional[List[str]] = None,
        custom_directives: Optional[List[str]] = None,
    ):
        """
        Initialize job script builder.
        
        Args:
            scheduler: Scheduler type (slurm, pbs, etc.)
            partition: Partition/queue name
            account: Account to charge
            time_limit: Job time limit
            modules: Environment modules to load
            custom_directives: Additional scheduler directives
        """
        self.scheduler = scheduler.lower()
        self.partition = partition
        self.account = account
        self.time_limit = time_limit
        self.modules = modules or []
        self.custom_directives = custom_directives or []
    
    def _get_partition_for_nodes(self, num_nodes: int, procs_per_node: int) -> str:
        """
        Select the best partition for a given number of nodes.
        
        Uses SlurmDetector to auto-select partition based on node count.
        Falls back to configured partition if detector unavailable.
        
        Args:
            num_nodes: Number of nodes for the job
            procs_per_node: Number of processes per node
            
        Returns:
            Partition name to use for the job
        """
        # If no detector available, use configured partition
        if not HAS_SLURM_DETECTOR or not self.partition:
            return self.partition or "X_usr_prod"
        
        try:
            detector = get_slurm_detector()
            if detector:
                best_partition = detector.select_best_partition(
                    num_nodes=num_nodes,
                    num_cpus_per_node=procs_per_node,
                    walltime_minutes=int(self.time_limit.split(":")[0]) * 60 + int(self.time_limit.split(":")[1])
                )
                if best_partition:
                    logger.debug(f"Auto-selected partition '{best_partition}' for {num_nodes} nodes")
                    return best_partition
        except Exception as e:
            logger.debug(f"Partition auto-selection failed: {e}")
        
        # Fallback to configured partition
        return self.partition or "X_usr_prod"
    
    def build_slurm_script(
        self,
        job_config: JobConfig,
        job_name: str,
        output_dir: Path,
        executable: Path,
        input_file: Path,
        launcher: str = "srun",
    ) -> str:
        """
        Build SLURM job script.
        
        Args:
            job_config: Job configuration
            job_name: Name of the job
            output_dir: Directory for output files
            executable: Path to executable
            input_file: Path to input file
            launcher: MPI launcher (srun, mpirun, mpiexec)
            
        Returns:
            Job script content as string
        """
        # Calculate procs per node
        procs_per_node = job_config.total_procs() // job_config.num_nodes if job_config.num_nodes > 0 else 1
        
        # Determine partition (auto-select if available)
        partition = self._get_partition_for_nodes(job_config.num_nodes, procs_per_node)
        
        script_lines = [
            "#!/bin/bash",
            "",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --nodes={job_config.num_nodes}",
            f"#SBATCH --ntasks={job_config.total_procs()}",
            f"#SBATCH --time={self.time_limit}",
            f"#SBATCH --output={output_dir / 'slurm-%j.out'}",
            f"#SBATCH --error={output_dir / 'slurm-%j.err'}",
        ]
        
        if partition:
            script_lines.append(f"#SBATCH --partition={partition}")
        
        if self.account:
            script_lines.append(f"#SBATCH --account={self.account}")
        
        # Add custom directives
        for directive in self.custom_directives:
            script_lines.append(f"#SBATCH {directive}")
        
        script_lines.extend([
            "",
            "# Environment setup",
            "module purge",
        ])
        
        # Load modules
        for module in self.modules:
            script_lines.append(f"module load {module}")
        
        script_lines.extend([
            "",
            "# Job information",
            "echo 'Job started at: $(date)'",
            f"echo 'Running on nodes: {job_config.num_nodes}'",
            f"echo 'Total processes: {job_config.total_procs()}'",
            f"echo 'MPI decomposition: {job_config.procs_decomposition}'",
            "",
            "# Run application",
            f"cd {output_dir}",
            f"{launcher} -n {job_config.total_procs()} {executable} {input_file.name}",
            "",
            "echo 'Job completed at: $(date)'",
        ])
        
        return "\n".join(script_lines)
    
    def build_pbs_script(
        self,
        job_config: JobConfig,
        job_name: str,
        output_dir: Path,
        executable: Path,
        input_file: Path,
        launcher: str = "mpirun",
    ) -> str:
        """
        Build PBS job script.
        
        Args:
            job_config: Job configuration
            job_name: Name of the job
            output_dir: Directory for output files
            executable: Path to executable
            input_file: Path to input file
            launcher: MPI launcher
            
        Returns:
            Job script content as string
        """
        script_lines = [
            "#!/bin/bash",
            "",
            f"#PBS -N {job_name}",
            f"#PBS -l nodes={job_config.num_nodes}",
            f"#PBS -l walltime={self.time_limit}",
            f"#PBS -o {output_dir / 'pbs.out'}",
            f"#PBS -e {output_dir / 'pbs.err'}",
        ]
        
        if self.partition:
            script_lines.append(f"#PBS -q {self.partition}")
        
        if self.account:
            script_lines.append(f"#PBS -A {self.account}")
        
        # Add custom directives
        for directive in self.custom_directives:
            script_lines.append(f"#PBS {directive}")
        
        script_lines.extend([
            "",
            "# Environment setup",
            "module purge",
        ])
        
        # Load modules
        for module in self.modules:
            script_lines.append(f"module load {module}")
        
        script_lines.extend([
            "",
            "# Change to working directory",
            "cd $PBS_O_WORKDIR",
            "",
            "# Job information",
            "echo 'Job started at: $(date)'",
            f"echo 'Running on nodes: {job_config.num_nodes}'",
            f"echo 'Total processes: {job_config.total_procs()}'",
            "",
            "# Run application",
            f"{launcher} -np {job_config.total_procs()} {executable} {input_file.name}",
            "",
            "echo 'Job completed at: $(date)'",
        ])
        
        return "\n".join(script_lines)
    
    def build(
        self,
        job_config: JobConfig,
        job_name: str,
        output_dir: Path,
        executable: Path,
        input_file: Path,
        launcher: str = "srun",
    ) -> str:
        """
        Build job script based on scheduler type.
        
        Args:
            job_config: Job configuration
            job_name: Name of the job
            output_dir: Directory for output files
            executable: Path to executable
            input_file: Path to input file
            launcher: MPI launcher
            
        Returns:
            Job script content as string
            
        Raises:
            ValueError: If scheduler type is not supported
        """
        if self.scheduler == "slurm":
            return self.build_slurm_script(
                job_config, job_name, output_dir, executable, input_file, launcher
            )
        elif self.scheduler == "pbs":
            return self.build_pbs_script(
                job_config, job_name, output_dir, executable, input_file, launcher
            )
        else:
            raise ValueError(f"Unsupported scheduler: {self.scheduler}")
    
    def write_script(
        self,
        script_content: str,
        output_path: Path,
    ) -> None:
        """
        Write job script to file.
        
        Args:
            script_content: Job script content
            output_path: Path to write script to
        """
        output_path.write_text(script_content)
        output_path.chmod(0o755)  # Make executable
        logger.info(f"Job script written to {output_path}")
