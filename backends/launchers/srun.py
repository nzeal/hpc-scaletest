"""
Slurm srun launcher backend with GPU binding support.
Optimized for Leonardo Booster and other GPU systems.
"""

import logging
from typing import List

from core.abstracts import LauncherInterface
from core.config import JobConfig, ResourceConfig


logger = logging.getLogger(__name__)


class SrunLauncher(LauncherInterface):
    """Slurm srun MPI launcher with GPU support."""
    
    def supports_gpu_binding(self) -> bool:
        """Check if launcher supports GPU binding."""
        return True  # srun has native GPU binding support via --gpu-bind
    
    def generate_launch_command(
        self,
        job_config: JobConfig,
        executable: List[str],
        resource_config: ResourceConfig
    ) -> List[str]:
        """Generate srun launch command with GPU binding."""
        cmd = ["srun"]
        
        # Resource specification
        # NOTE: Don't specify -n/--ntasks here - let SLURM use #SBATCH --ntasks
        # Only specify per-node allocation
        cmd.extend(["--ntasks-per-node", str(resource_config.procs_per_node)])
        
        # GPU binding - critical for GPU jobs
        if resource_config.gpus_per_node > 0:
            cmd.extend(["--gpus-per-node", str(resource_config.gpus_per_node)])
            cmd.append("--gpu-bind=closest")
            logger.debug(f"srun GPU config: {resource_config.gpus_per_node} GPUs/node, "
                       f"binding=closest, {resource_config.procs_per_node} tasks/node")
            
            # Verify 1:1 ratio
            if resource_config.procs_per_node != resource_config.gpus_per_node:
                logger.warning(f"⚠ GPU/task mismatch: {resource_config.gpus_per_node} GPUs but "
                             f"{resource_config.procs_per_node} tasks per node")
        
        # CPU binding
        cmd.append("--cpu-bind=cores")
        
        # Custom options
        for key, value in self.options.items():
            if value is True:
                cmd.append(f"--{key}")
            elif value is not False:
                cmd.extend([f"--{key}", str(value)])
        
        # Executable
        cmd.extend(executable)
        
        return cmd
