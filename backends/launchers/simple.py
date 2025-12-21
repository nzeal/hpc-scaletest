"""
Simple launcher backend that doesn't add MPI options.
"""

from typing import List

from core.abstracts import LauncherInterface
from core.config import JobConfig, ResourceConfig


class SimpleLauncher(LauncherInterface):
    """Simple launcher that just returns the executable as-is."""
    
    def generate_launch_command(
        self,
        job_config: JobConfig,
        executable: List[str],
        resource_config: ResourceConfig
    ) -> List[str]:
        """Return the executable without any launcher options."""
        return executable
    
    def supports_gpu_binding(self) -> bool:
        return False