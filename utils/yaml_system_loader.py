"""
YAML-based System Configuration Loader.

Loads system configurations from external YAML files, removing hardcoded
assumptions about specific HPC systems.
"""

import logging
import os
import socket
from pathlib import Path
from typing import Optional, Dict, Any, List

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False
    yaml = None

logger = logging.getLogger(__name__)


class SystemConfig:
    """Configuration for a specific HPC system."""
    
    def __init__(self, name: str, config: Dict[str, Any]):
        self.name = name
        self.description = config.get('description', '')
        self.scheduler = config.get('scheduler', 'slurm')
        self.launcher = config.get('launcher', 'srun')
        self.module_system = config.get('module_system', 'lmod')
        self.partitions = config.get('partitions', {})
        self.partition_aliases = config.get('partition_aliases', {})
    
    def get_partition(self, name: str) -> Optional[Dict[str, Any]]:
        """Get partition configuration, resolving aliases."""
        # Check for alias first
        if name in self.partition_aliases:
            name = self.partition_aliases[name]
        
        return self.partitions.get(name)
    
    def list_partitions(self) -> List[str]:
        """List all available partition names."""
        return list(self.partitions.keys())


class YAMLSystemLoader:
    """
    Load system configurations from YAML files.
    
    Searches for configuration files in standard locations and allows
    users to define their own HPC system configurations.
    """
    
    DEFAULT_CONFIG_PATHS = [
        Path(__file__).parent.parent / 'config' / 'systems.yaml',
        Path.home() / '.hpc-scaletest' / 'systems.yaml',
        Path('/etc/hpc-scaletest/systems.yaml'),
    ]
    
    def __init__(self, config_path: Optional[Path] = None):
        """
        Initialize the system loader.
        
        Args:
            config_path: Path to YAML configuration file.
                        If None, searches default locations.
        """
        self.config_path = config_path
        self.systems: Dict[str, SystemConfig] = {}
        self.mpi_implementations: Dict[str, Dict] = {}
        self.gpu_vendors: Dict[str, Dict] = {}
        self._loaded = False
        
        # Load configuration
        self._load_config()
    
    def _find_config_file(self) -> Optional[Path]:
        """Find configuration file in default locations."""
        if self.config_path and self.config_path.exists():
            return self.config_path
        
        # Check environment variable
        env_path = os.environ.get('HPC_SCALETEST_CONFIG')
        if env_path:
            path = Path(env_path)
            if path.exists():
                return path
        
        # Search default paths
        for path in self.DEFAULT_CONFIG_PATHS:
            if path.exists():
                logger.debug(f"Found config file: {path}")
                return path
        
        return None
    
    def _load_config(self) -> bool:
        """Load configuration from YAML file."""
        if not HAS_YAML:
            logger.warning("PyYAML not installed, using built-in defaults")
            self._load_builtin_defaults()
            return False
        
        config_file = self._find_config_file()
        if not config_file:
            logger.info("No external config file found, using built-in defaults")
            self._load_builtin_defaults()
            return False
        
        try:
            with open(config_file, 'r') as f:
                config = yaml.safe_load(f)
            
            # Load systems
            for name, system_config in config.get('systems', {}).items():
                self.systems[name] = SystemConfig(name, system_config)
            
            # Load MPI implementations
            self.mpi_implementations = config.get('mpi_implementations', {})
            
            # Load GPU vendors
            self.gpu_vendors = config.get('gpu_vendors', {})
            
            self._loaded = True
            logger.info(f"Loaded {len(self.systems)} system configurations from {config_file}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to load config from {config_file}: {e}")
            self._load_builtin_defaults()
            return False
    
    def _load_builtin_defaults(self):
        """Load minimal built-in defaults when no config file is available."""
        # Generic system that works with auto-detection
        self.systems['generic'] = SystemConfig('generic', {
            'description': 'Generic HPC system with auto-detection',
            'scheduler': 'slurm',
            'launcher': 'srun',
            'module_system': 'lmod',
            'partitions': {
                'default': {
                    'cores_per_node': 'auto',
                    'gpus_per_node': 0,
                    'memory_gb': 128,
                    'time_limit': '01:00:00'
                }
            }
        })
        
        # MPI defaults
        self.mpi_implementations = {
            'openmpi': {
                'supports_map_by': True,
                'supports_bind_to': True,
                'supports_report_bindings': True,
                'launcher': 'mpirun'
            },
            'intelmpi': {
                'supports_map_by': True,
                'supports_bind_to': True,
                'supports_report_bindings': False,
                'launcher': 'mpirun'
            },
            'mpich': {
                'supports_map_by': False,
                'supports_bind_to': False,
                'supports_report_bindings': False,
                'launcher': 'mpiexec'
            }
        }
        
        # GPU vendor defaults
        self.gpu_vendors = {
            'nvidia': {'device_env_var': 'CUDA_VISIBLE_DEVICES'},
            'amd': {'device_env_var': 'ROCR_VISIBLE_DEVICES'},
            'intel': {'device_env_var': 'ZE_AFFINITY_MASK'}
        }
    
    def get_system(self, name: str) -> Optional[SystemConfig]:
        """Get a system configuration by name."""
        return self.systems.get(name)
    
    def detect_system(self) -> Optional[SystemConfig]:
        """
        Auto-detect the current system based on hostname or environment.
        
        Returns:
            SystemConfig if detected, None otherwise
        """
        # Check environment variable override
        env_system = os.environ.get('HPC_SCALETEST_SYSTEM')
        if env_system and env_system in self.systems:
            logger.info(f"Using system from environment: {env_system}")
            return self.systems[env_system]
        
        # Try to detect from hostname
        hostname = socket.gethostname().lower()
        
        # Simple hostname matching
        hostname_hints = {
            'leonardo': 'leonardo_booster',
            'lumi': 'lumi',
            'dgx': 'dgx_a100',
        }
        
        for hint, system_name in hostname_hints.items():
            if hint in hostname and system_name in self.systems:
                logger.info(f"Auto-detected system: {system_name}")
                return self.systems[system_name]
        
        # Fall back to generic
        logger.info("Using generic system configuration")
        return self.systems.get('generic')
    
    def list_systems(self) -> List[str]:
        """List all available system names."""
        return list(self.systems.keys())
    
    def get_mpi_config(self, implementation: str) -> Dict[str, Any]:
        """Get MPI implementation configuration."""
        return self.mpi_implementations.get(implementation, {})
    
    def get_gpu_vendor_config(self, vendor: str) -> Dict[str, Any]:
        """Get GPU vendor configuration."""
        return self.gpu_vendors.get(vendor.lower(), {})
    
    def get_gpu_device_env_var(self, vendor: str) -> str:
        """Get the environment variable name for GPU device selection."""
        config = self.get_gpu_vendor_config(vendor)
        return config.get('device_env_var', 'CUDA_VISIBLE_DEVICES')


# Global loader instance (lazy initialization)
_global_loader: Optional[YAMLSystemLoader] = None


def get_system_loader() -> YAMLSystemLoader:
    """Get the global system loader instance."""
    global _global_loader
    if _global_loader is None:
        _global_loader = YAMLSystemLoader()
    return _global_loader


def get_partition_config(system_name: str, partition_name: str) -> Optional[Dict[str, Any]]:
    """
    Convenience function to get partition configuration.
    
    Args:
        system_name: Name of the system (e.g., 'leonardo_booster')
        partition_name: Name of the partition (e.g., 'boost_usr_prod')
    
    Returns:
        Partition configuration dict or None
    """
    loader = get_system_loader()
    system = loader.get_system(system_name)
    if system:
        return system.get_partition(partition_name)
    return None
