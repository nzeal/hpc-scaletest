"""
Launcher registry system for custom launcher definitions.
"""

import logging
from typing import Dict, Type, Optional, List
from core.abstracts import LauncherInterface
from core.config import JobConfig, ResourceConfig


logger = logging.getLogger(__name__)


# Global launcher registry
_LAUNCHER_REGISTRY: Dict[str, Type[LauncherInterface]] = {}


def register_launcher(name: str):
    """
    Decorator to register a custom launcher.
    
    Example:
        @register_launcher('mpirun-mapby')
        class MpirunMapbyLauncher(LauncherInterface):
            def command(self, job):
                return ['mpirun', '-np', str(job.num_procs)]
    """
    def decorator(cls):
        if not issubclass(cls, LauncherInterface):
            raise TypeError(f"{cls.__name__} must inherit from LauncherInterface")
        
        if name in _LAUNCHER_REGISTRY:
            logger.warning(f"Overwriting existing launcher registration: {name}")
        
        _LAUNCHER_REGISTRY[name] = cls
        logger.debug(f"Registered launcher: {name} -> {cls.__name__}")
        return cls
    
    return decorator


def get_launcher(name: str, options: Optional[Dict] = None) -> LauncherInterface:
    """
    Get a launcher instance by name from the registry.
    
    Args:
        name: Launcher name
        options: Optional configuration options
        
    Returns:
        LauncherInterface instance
        
    Raises:
        ValueError: If launcher not found
    """
    if name not in _LAUNCHER_REGISTRY:
        raise ValueError(
            f"Launcher '{name}' not found in registry. "
            f"Available launchers: {list(_LAUNCHER_REGISTRY.keys())}"
        )
    
    launcher_cls = _LAUNCHER_REGISTRY[name]
    return launcher_cls(options)


def list_launchers() -> List[str]:
    """List all registered launcher names."""
    return list(_LAUNCHER_REGISTRY.keys())


def has_launcher(name: str) -> bool:
    """Check if a launcher is registered."""
    return name in _LAUNCHER_REGISTRY


class JobLauncher(LauncherInterface):
    """
    Base class for custom job launchers.
    Simplifies launcher creation by providing a command() method interface.
    """
    
    def command(self, job: JobConfig, resource_config: ResourceConfig) -> List[str]:
        """
        Generate launch command for the job.
        
        Args:
            job: Job configuration
            resource_config: Resource configuration
            
        Returns:
            List of command arguments
        """
        raise NotImplementedError("Subclasses must implement command()")
    
    def generate_launch_command(
        self,
        job_config: JobConfig,
        executable: List[str],
        resource_config: ResourceConfig
    ) -> List[str]:
        """Generate the MPI launch command."""
        # Get the launcher prefix from command()
        cmd = self.command(job_config, resource_config)
        # Append the executable and its arguments
        cmd.extend(executable)
        return cmd
    
    def supports_gpu_binding(self) -> bool:
        """Check if launcher supports GPU binding."""
        return self.options.get('supports_gpu_binding', False)
