"""
Slurm srun launcher backend with GPU binding support.
"""

from typing import List

from core.abstracts import LauncherInterface
from core.config import JobConfig, ResourceConfig


class SrunLauncher(LauncherInterface):
    """Slurm srun MPI launcher."""
    
    def generate_launch_command(
        self,
        job_config: JobConfig,
        executable: List[str],
        resource_config: ResourceConfig
    ) -> List[str]:
        """Generate srun launch command."""
        cmd = ["srun"]
        
        # Resource specification
        cmd.extend(["--ntasks-per-node", str(resource_config.procs_per_node)])
        
        # GPU binding
        if resource_config.gpus_per_node > 0:
            cmd.extend(["--gpus-per-node", str(resource_config.gpus_per_node)])
            cmd.append("--gpu-bind=closest")
        
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
    
    def supports_gpu_binding(self) -> bool:
        return True
