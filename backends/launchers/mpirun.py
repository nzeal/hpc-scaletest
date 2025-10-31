"""
Generic mpirun/mpiexec launcher backend.
"""

from typing import List

from core.abstracts import LauncherInterface
from core.config import JobConfig, ResourceConfig


class MpiRunLauncher(LauncherInterface):
    """Generic MPI launcher (mpirun/mpiexec)."""
    
    def generate_launch_command(
        self,
        job_config: JobConfig,
        executable: List[str],
        resource_config: ResourceConfig
    ) -> List[str]:
        """Generate mpirun launch command."""
        launcher = self.options.get('launcher', 'mpirun')
        cmd = [launcher, "-np", str(job_config.num_procs)]
        
        # Procs per node
        if resource_config.procs_per_node:
            cmd.extend(["--npernode", str(resource_config.procs_per_node)])
        
        # Mapping/binding (default for template)
        map_by = self.options.get('map_by', 'ppr:4:node:PE=8')
        cmd.extend(["--map-by", map_by])
        if self.options.get('report_bindings', True):
            cmd.append("--report-bindings")
        
        # Custom options
        for key, value in self.options.items():
            if key not in ('launcher', 'map_by', 'report_bindings'):
                if value is True:
                    cmd.append(f"--{key}")
                elif value is not False:
                    cmd.extend([f"--{key}", str(value)])
        
        # Executable
        cmd.extend(executable)
        
        return cmd
    
    def supports_gpu_binding(self) -> bool:
        return False
