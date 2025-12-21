"""
Core decorators for framework extensibility.

Provides decorators for registering:
- Scheduler directives
- Hardware features
- Scheduler backends
- MPI launchers
- Module systems
- Custom validators

All decorators use the central PluginRegistry for storage.
"""

import logging
from functools import wraps
from core.registry import PluginRegistry

logger = logging.getLogger(__name__)


def register_directive(name):
    """
    Decorator to register a scheduler directive function.
    
    Example:
        @register_directive("gpu")
        def gpu_directive(num_gpus, **kwargs):
            return {'slurm': '--gpus-per-node=4', ...}
    
    Args:
        name: Directive name (e.g., "gpu", "memory", "network")
    """
    def decorator(func):
        PluginRegistry.register_directive(name, func)
        
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        
        return wrapper
    
    return decorator


def register_feature(name):
    """
    Decorator to register a hardware feature class.
    
    Example:
        @register_feature("hbm")
        class HighBandwidthMemory(HardwareFeature):
            def detect(self):
                return True
    
    Args:
        name: Feature name (e.g., "hbm", "nvme", "infiniband")
    """
    def decorator(cls):
        PluginRegistry.register_feature(name, cls)
        return cls
    
    return decorator


def register_backend(name):
    """
    Decorator to register a scheduler backend class.
    
    Example:
        @register_backend("slurm")
        class SlurmScheduler(SchedulerInterface):
            ...
    
    Args:
        name: Backend name (e.g., "slurm", "pbs", "lsf")
    """
    def decorator(cls):
        PluginRegistry.register_backend(name, cls)
        return cls
    
    return decorator


def register_module_system(name):
    """
    Decorator to register a module system class.
    
    Example:
        @register_module_system("lmod")
        class LmodSystem(ModuleSystemInterface):
            ...
    
    Args:
        name: Module system name (e.g., "lmod", "environment-modules")
    """
    def decorator(cls):
        PluginRegistry.register_module_system(name, cls)
        return cls
    
    return decorator


def register_validator(name):
    """
    Decorator to register a configuration validator.
    
    Example:
        @register_validator("gpu_config")
        def validate_gpu_config(config):
            return True, []
    
    Args:
        name: Validator name
    """
    def decorator(func):
        PluginRegistry.register_validator(name, func)
        
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        
        return wrapper
    
    return decorator


# Convenience getters - delegate to PluginRegistry

def get_directive(name):
    """Get registered directive function by name."""
    directive = PluginRegistry.get_directive(name)
    if directive is None:
        raise KeyError(f"Directive '{name}' not registered")
    return directive


def get_feature(name):
    """Get registered feature class by name."""
    feature = PluginRegistry.get_feature(name)
    if feature is None:
        raise KeyError(f"Feature '{name}' not registered")
    return feature


def get_backend(name):
    """Get registered backend class by name."""
    backend = PluginRegistry.get_backend(name)
    if backend is None:
        raise KeyError(f"Backend '{name}' not registered")
    return backend


def get_validator(name):
    """Get registered validator function by name."""
    validator = PluginRegistry.get_validator(name)
    if validator is None:
        raise KeyError(f"Validator '{name}' not registered")
    return validator


def list_directives():
    """List all registered directive names."""
    return PluginRegistry.list_directives()


def list_features():
    """List all registered feature names."""
    return PluginRegistry.list_features()


def list_backends():
    """List all registered backend names."""
    return PluginRegistry.list_backends()


def list_validators():
    """List all registered validator names."""
    return PluginRegistry.list_validators()

