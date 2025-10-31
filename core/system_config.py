"""
System configuration structures for HPC systems.

Similar to ReFrame's system configuration, this allows defining:
- System information (hostnames, partitions)
- Resources per partition (CPUs, GPUs, memory)
- Environments (modules, compilers)
- Logging and storage backends
"""

import json
import re
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Dict, Optional, Any
from importlib.util import spec_from_file_location, module_from_spec


logger = logging.getLogger(__name__)


@dataclass
class ProcessorConfig:
    """Processor configuration for a partition."""
    arch: str = "unknown"
    platform: str = "x86_64"
    num_cpus: int = 1
    num_cpus_per_core: int = 1
    num_cpus_per_socket: int = 1
    num_sockets: int = 1


@dataclass
class DeviceConfig:
    """Device configuration (e.g., GPUs)."""
    type: str = "gpu"
    arch: str = "unknown"
    model: str = "unknown"
    num_devices: int = 0


@dataclass
class ContainerPlatformConfig:
    """Container platform configuration."""
    type: str = "Singularity"
    default: bool = True
    modules: List[str] = field(default_factory=list)
    env_vars: List[List[str]] = field(default_factory=list)


@dataclass
class PartitionConfig:
    """Partition configuration within a system."""
    name: str
    descr: str
    scheduler: str
    launcher: str
    modules: List[str] = field(default_factory=list)
    access: List[str] = field(default_factory=list)
    resources: List[Dict[str, Any]] = field(default_factory=list)
    max_jobs: int = 1
    environs: List[str] = field(default_factory=list)
    processor: Optional[Dict[str, Any]] = None
    devices: List[Dict[str, Any]] = field(default_factory=list)
    extras: Dict[str, Any] = field(default_factory=dict)
    container_platforms: List[Dict[str, Any]] = field(default_factory=list)
    sched_options: Dict[str, Any] = field(default_factory=dict)


@dataclass
class SystemConfig:
    """System configuration."""
    name: str
    descr: str
    hostnames: List[str] = field(default_factory=list)
    modules_system: str = "nomod"
    partitions: List[PartitionConfig] = field(default_factory=list)


@dataclass
class EnvironmentDef:
    """Environment definition with compilers and modules."""
    name: str
    modules: List[str] = field(default_factory=list)
    cc: str = "gcc"
    cxx: str = "g++"
    ftn: str = "gfortran"
    nvcc: Optional[str] = None
    features: List[str] = field(default_factory=list)
    extras: Dict[str, Any] = field(default_factory=dict)


@dataclass
class LoggingHandlerConfig:
    """Logging handler configuration."""
    type: str
    name: Optional[str] = None
    level: str = "info"
    format: str = "%(message)s"
    append: bool = False
    prefix: Optional[str] = None
    format_perfvars: Optional[str] = None


@dataclass
class LoggingConfig:
    """Logging configuration."""
    level: str = "info"
    handlers: List[Dict[str, Any]] = field(default_factory=list)
    handlers_perflog: List[Dict[str, Any]] = field(default_factory=list)


@dataclass
class StorageConfig:
    """Storage backend configuration."""
    enable: bool = False
    backend: str = "postgresql"
    postgresql_driver: Optional[str] = None
    postgresql_host: Optional[str] = None
    postgresql_port: Optional[int] = None
    postgresql_db: Optional[str] = None
    postgresql_conn_timeout: int = 60


@dataclass
class ModeConfig:
    """Execution mode configuration."""
    name: str
    options: List[str] = field(default_factory=list)


@dataclass
class GeneralConfig:
    """General configuration."""
    check_search_path: List[str] = field(default_factory=lambda: ["checks/"])
    check_search_recursive: bool = True


@dataclass
class SiteConfiguration:
    """Complete site configuration."""
    systems: List[SystemConfig] = field(default_factory=list)
    environments: List[EnvironmentDef] = field(default_factory=list)
    logging: List[LoggingConfig] = field(default_factory=list)
    storage: List[StorageConfig] = field(default_factory=list)
    modes: List[ModeConfig] = field(default_factory=list)
    general: List[GeneralConfig] = field(default_factory=list)
    
    def get_system(self, name: str) -> Optional[SystemConfig]:
        """Get system by name."""
        for system in self.systems:
            if system.name == name:
                return system
        return None
    
    def get_environment(self, name: str) -> Optional[EnvironmentDef]:
        """Get environment by name."""
        for env in self.environments:
            if env.name == name:
                return env
        return None
    
    def find_system_by_hostname(self, hostname: str) -> Optional[SystemConfig]:
        """Find system matching the given hostname."""
        for system in self.systems:
            for pattern in system.hostnames:
                if re.match(pattern, hostname):
                    return system
        return None


def load_system_config(config_file: Path) -> Optional[SiteConfiguration]:
    """
    Load system configuration from a Python file.
    
    The file should define a `site_configuration` dictionary with the structure
    similar to ReFrame's configuration.
    
    Args:
        config_file: Path to the system configuration Python file
        
    Returns:
        SiteConfiguration object or None if loading fails
    """
    try:
        # Import the configuration module
        spec = spec_from_file_location("system_config", config_file)
        if spec is None or spec.loader is None:
            logger.error(f"Failed to load spec for {config_file}")
            return None
        
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
        
        # Get the site_configuration dictionary
        if not hasattr(module, 'site_configuration'):
            logger.error(f"No 'site_configuration' found in {config_file}")
            return None
        
        config_dict = module.site_configuration
        
        # Parse the configuration
        return parse_site_configuration(config_dict)
        
    except Exception as e:
        logger.error(f"Failed to load system configuration from {config_file}: {e}")
        return None


def parse_site_configuration(config_dict: Dict) -> SiteConfiguration:
    """Parse site configuration dictionary into dataclass structure."""
    
    site_config = SiteConfiguration()
    
    # Parse systems
    for system_dict in config_dict.get('systems', []):
        partitions = []
        for part_dict in system_dict.get('partitions', []):
            partition = PartitionConfig(
                name=part_dict['name'],
                descr=part_dict.get('descr', ''),
                scheduler=part_dict.get('scheduler', 'local'),
                launcher=part_dict.get('launcher', 'local'),
                modules=part_dict.get('modules', []),
                access=part_dict.get('access', []),
                resources=part_dict.get('resources', []),
                max_jobs=part_dict.get('max_jobs', 1),
                environs=part_dict.get('environs', []),
                processor=part_dict.get('processor'),
                devices=part_dict.get('devices', []),
                extras=part_dict.get('extras', {}),
                container_platforms=part_dict.get('container_platforms', []),
                sched_options=part_dict.get('sched_options', {})
            )
            partitions.append(partition)
        
        system = SystemConfig(
            name=system_dict['name'],
            descr=system_dict.get('descr', ''),
            hostnames=system_dict.get('hostnames', []),
            modules_system=system_dict.get('modules_system', 'nomod'),
            partitions=partitions
        )
        site_config.systems.append(system)
    
    # Parse environments
    for env_dict in config_dict.get('environments', []):
        env = EnvironmentDef(
            name=env_dict['name'],
            modules=env_dict.get('modules', []),
            cc=env_dict.get('cc', 'gcc'),
            cxx=env_dict.get('cxx', 'g++'),
            ftn=env_dict.get('ftn', 'gfortran'),
            nvcc=env_dict.get('nvcc'),
            features=env_dict.get('features', []),
            extras=env_dict.get('extras', {})
        )
        site_config.environments.append(env)
    
    # Parse logging
    for log_dict in config_dict.get('logging', []):
        log_config = LoggingConfig(
            level=log_dict.get('level', 'info'),
            handlers=log_dict.get('handlers', []),
            handlers_perflog=log_dict.get('handlers_perflog', [])
        )
        site_config.logging.append(log_config)
    
    # Parse storage
    for storage_dict in config_dict.get('storage', []):
        storage = StorageConfig(
            enable=storage_dict.get('enable', False),
            backend=storage_dict.get('backend', 'postgresql'),
            postgresql_driver=storage_dict.get('postgresql_driver'),
            postgresql_host=storage_dict.get('postgresql_host'),
            postgresql_port=storage_dict.get('postgresql_port'),
            postgresql_db=storage_dict.get('postgresql_db'),
            postgresql_conn_timeout=storage_dict.get('postgresql_conn_timeout', 60)
        )
        site_config.storage.append(storage)
    
    # Parse modes
    for mode_dict in config_dict.get('modes', []):
        mode = ModeConfig(
            name=mode_dict['name'],
            options=mode_dict.get('options', [])
        )
        site_config.modes.append(mode)
    
    # Parse general
    for general_dict in config_dict.get('general', []):
        general = GeneralConfig(
            check_search_path=general_dict.get('check_search_path', ['checks/']),
            check_search_recursive=general_dict.get('check_search_recursive', True)
        )
        site_config.general.append(general)
    
    return site_config


def format_record(record: Any, extras: Dict, ignore_keys: List[str]) -> str:
    """
    Format a report record to match the structure expected by the database.
    
    Custom JSON formatter that splits pipe-separated check_perf_var
    into individual fields before sending to HTTP endpoint.
    """
    
    def _sanitize(s):
        return re.sub(r'\W', '_', s)
    
    json_record = {}
    
    formatted_record = {
        'session_info': {
            'uuid': record.get('session_info', {}).get('uuid'),
            'time_start_unix': record.get('session_info', {}).get('time_start_unix'),
            'time_end_unix': record.get('session_info', {}).get('time_end_unix'),
            'time_elapsed': (
                record.get('session_info', {}).get('time_end_unix', 0) - 
                record.get('session_info', {}).get('time_start_unix', 0)
            ),
        },
        'runs': []
    }

    # Handle dict or object with __dict__
    if hasattr(record, '__dict__'):
        record_dict = record.__dict__
    else:
        record_dict = record

    for attr, val in record_dict.items():
        if attr in ignore_keys or attr.startswith('_'):
            continue

        # Special handling to convert MappingView or other non-serializable types
        if attr == 'check_perfvalues':
            # Convert to plain dict if needed
            val = dict(val)
            sanitized_val = {}
            for k, v in val.items():
                _sanitized_key = _sanitize(k.split(':')[-1])
                sanitized_val[_sanitized_key] = v
            val = sanitized_val

        if attr.startswith('check_'):
            json_record[attr[6:]] = val
        else:
            json_record[attr] = val

    json_record.update(extras)
    
    return json.dumps(formatted_record)
