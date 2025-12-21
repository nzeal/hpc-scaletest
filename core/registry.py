"""
Plugin Registry System

Central registry for all plugin types in HPC-ScaleTest.
Enables runtime discovery and loading of plugins without
modifying core code.

Plugin Types:
    - Scheduler Directives
    - Hardware Features
    - Scheduler Backends
    - MPI Launchers
    - Module Systems
    - Validators
"""

import logging
from typing import Dict, Type, Optional, List, Callable, Any
from core.abstracts import LauncherInterface
from core.config import JobConfig, ResourceConfig


logger = logging.getLogger(__name__)


class PluginRegistry:
    """
    Central registry for all plugin types.
    
    Uses separate registries for each plugin type to avoid
    naming conflicts and provide type-specific operations.
    """
    
    # Plugin type registries
    _directives: Dict[str, Callable] = {}
    _features: Dict[str, Type] = {}
    _backends: Dict[str, Type] = {}
    _launchers: Dict[str, Type[LauncherInterface]] = {}
    _module_systems: Dict[str, Type] = {}
    _validators: Dict[str, Callable] = {}
    
    # Metadata storage
    _metadata: Dict[str, Dict[str, Any]] = {}
    
    @classmethod
    def register_directive(cls, name: str, func: Callable, metadata: Optional[Dict] = None):
        """Register a scheduler directive function."""
        if name in cls._directives:
            logger.warning(f"Directive '{name}' already registered, overwriting")
        cls._directives[name] = func
        if metadata:
            cls._metadata[f"directive:{name}"] = metadata
        logger.debug(f"Registered directive: {name}")
    
    @classmethod
    def register_feature(cls, name: str, feature_class: Type, metadata: Optional[Dict] = None):
        """Register a hardware feature class."""
        if name in cls._features:
            logger.warning(f"Feature '{name}' already registered, overwriting")
        cls._features[name] = feature_class
        if metadata:
            cls._metadata[f"feature:{name}"] = metadata
        logger.debug(f"Registered feature: {name}")
    
    @classmethod
    def register_backend(cls, name: str, backend_class: Type, metadata: Optional[Dict] = None):
        """Register a scheduler backend class."""
        if name in cls._backends:
            logger.warning(f"Backend '{name}' already registered, overwriting")
        cls._backends[name] = backend_class
        if metadata:
            cls._metadata[f"backend:{name}"] = metadata
        logger.debug(f"Registered backend: {name}")
    
    @classmethod
    def register_launcher(cls, name: str, launcher_class: Type[LauncherInterface], metadata: Optional[Dict] = None):
        """Register an MPI launcher class."""
        if not issubclass(launcher_class, LauncherInterface):
            raise TypeError(f"{launcher_class.__name__} must inherit from LauncherInterface")
        if name in cls._launchers:
            logger.warning(f"Launcher '{name}' already registered, overwriting")
        cls._launchers[name] = launcher_class
        if metadata:
            cls._metadata[f"launcher:{name}"] = metadata
        logger.debug(f"Registered launcher: {name}")
    
    @classmethod
    def register_module_system(cls, name: str, module_class: Type, metadata: Optional[Dict] = None):
        """Register a module system class."""
        if name in cls._module_systems:
            logger.warning(f"Module system '{name}' already registered, overwriting")
        cls._module_systems[name] = module_class
        if metadata:
            cls._metadata[f"module:{name}"] = metadata
        logger.debug(f"Registered module system: {name}")
    
    @classmethod
    def register_validator(cls, name: str, func: Callable, metadata: Optional[Dict] = None):
        """Register a validator function."""
        if name in cls._validators:
            logger.warning(f"Validator '{name}' already registered, overwriting")
        cls._validators[name] = func
        if metadata:
            cls._metadata[f"validator:{name}"] = metadata
        logger.debug(f"Registered validator: {name}")
    
    # Getter methods
    @classmethod
    def get_directive(cls, name: str) -> Optional[Callable]:
        return cls._directives.get(name)
    
    @classmethod
    def get_feature(cls, name: str) -> Optional[Type]:
        return cls._features.get(name)
    
    @classmethod
    def get_backend(cls, name: str) -> Optional[Type]:
        return cls._backends.get(name)
    
    @classmethod
    def get_launcher(cls, name: str, options: Optional[Dict] = None) -> Optional[LauncherInterface]:
        """Get launcher instance by name."""
        launcher_cls = cls._launchers.get(name)
        if launcher_cls is None:
            return None
        return launcher_cls(options)
    
    @classmethod
    def get_module_system(cls, name: str) -> Optional[Type]:
        return cls._module_systems.get(name)
    
    @classmethod
    def get_validator(cls, name: str) -> Optional[Callable]:
        return cls._validators.get(name)
    
    # List methods
    @classmethod
    def list_directives(cls) -> List[str]:
        return list(cls._directives.keys())
    
    @classmethod
    def list_features(cls) -> List[str]:
        return list(cls._features.keys())
    
    @classmethod
    def list_backends(cls) -> List[str]:
        return list(cls._backends.keys())
    
    @classmethod
    def list_launchers(cls) -> List[str]:
        return list(cls._launchers.keys())
    
    @classmethod
    def list_module_systems(cls) -> List[str]:
        return list(cls._module_systems.keys())
    
    @classmethod
    def list_validators(cls) -> List[str]:
        return list(cls._validators.keys())
    
    @classmethod
    def list_all_plugins(cls) -> Dict[str, List[str]]:
        """List all registered plugins by type."""
        return {
            'directives': cls.list_directives(),
            'features': cls.list_features(),
            'backends': cls.list_backends(),
            'launchers': cls.list_launchers(),
            'module_systems': cls.list_module_systems(),
            'validators': cls.list_validators(),
        }
    
    @classmethod
    def print_registry_status(cls):
        """Print current registry status (for debugging)."""
        print("="*70)
        print("HPC-ScaleTest Plugin Registry Status")
        print("="*70)
        plugins = cls.list_all_plugins()
        for plugin_type, names in plugins.items():
            count = len(names)
            print(f"\n{plugin_type.title()}: {count}")
            if names:
                for name in sorted(names):
                    print(f"  • {name}")
        print()
        print(f"Total plugins: {sum(len(names) for names in plugins.values())}")
        print("="*70)


# Legacy compatibility
_LAUNCHER_REGISTRY: Dict[str, Type[LauncherInterface]] = PluginRegistry._launchers

def register_launcher(name: str):
    """
    Decorator to register a custom launcher.
    
    Example:
        @register_launcher('mpirun-mapby')
        class MpirunMapbyLauncher(LauncherInterface):
            def generate_launch_command(self, job_config, executable, resource_config):
                return ['mpirun', '-np', str(resource_config.num_procs)]
    """
    def decorator(cls):
        PluginRegistry.register_launcher(name, cls)
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
    launcher = PluginRegistry.get_launcher(name, options)
    if launcher is None:
        raise ValueError(
            f"Launcher '{name}' not found in registry. "
            f"Available launchers: {PluginRegistry.list_launchers()}"
        )
    return launcher


def list_launchers() -> List[str]:
    """List all registered launcher names."""
    return PluginRegistry.list_launchers()


def has_launcher(name: str) -> bool:
    """Check if a launcher is registered."""
    return name in PluginRegistry._launchers



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
