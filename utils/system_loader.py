"""
System configuration loader utility.

This module provides utilities to load and apply system configurations
from Python configuration files.
"""

import logging
import socket
from pathlib import Path
from typing import Optional

from core.system_config import (
    SiteConfiguration, 
    SystemConfig, 
    PartitionConfig,
    load_system_config
)
from core.config import ResourceConfig, BackendConfig
from core.types import SchedulerBackend, LauncherBackend, ModuleBackend


logger = logging.getLogger(__name__)


class SystemConfigLoader:
    """
    Loader for system configurations.
    
    Automatically detects the current system based on hostname
    and loads the appropriate configuration.
    """
    
    def __init__(self, config_file: Optional[Path] = None):
        """
        Initialize the system configuration loader.
        
        Args:
            config_file: Path to the system configuration file.
                        If None, will search for common names like leonardo_system.py
        """
        self.config_file = config_file
        self.site_config: Optional[SiteConfiguration] = None
        self.current_system: Optional[SystemConfig] = None
        
        if config_file and config_file.exists():
            self.load_configuration(config_file)
    
    def load_configuration(self, config_file: Path) -> bool:
        """
        Load system configuration from file.
        
        Args:
            config_file: Path to the configuration file
            
        Returns:
            True if loaded successfully, False otherwise
        """
        logger.info(f"Loading system configuration from {config_file}")
        self.site_config = load_system_config(config_file)
        
        if self.site_config:
            # Try to detect current system
            hostname = socket.gethostname()
            self.current_system = self.site_config.find_system_by_hostname(hostname)
            
            if self.current_system:
                logger.info(f"Detected system: {self.current_system.name}")
            else:
                logger.warning(f"Hostname '{hostname}' does not match any system in configuration")
            
            return True
        
        return False
    
    def find_config_file(self, search_paths: list[Path]) -> Optional[Path]:
        """
        Search for a system configuration file in the given paths.
        
        Args:
            search_paths: List of directories to search
            
        Returns:
            Path to the configuration file if found, None otherwise
        """
        common_names = [
            'leonardo_system.py',
            'system_config.py',
            'site_config.py',
        ]
        
        for search_path in search_paths:
            for name in common_names:
                config_file = search_path / name
                if config_file.exists():
                    logger.info(f"Found system configuration: {config_file}")
                    return config_file
        
        return None
    
    def get_partition(self, partition_name: str) -> Optional[PartitionConfig]:
        """
        Get a partition configuration by name.
        
        Args:
            partition_name: Name of the partition
            
        Returns:
            PartitionConfig if found, None otherwise
        """
        if not self.current_system:
            return None
        
        for partition in self.current_system.partitions:
            if partition.name == partition_name:
                return partition
        
        return None
    
    def create_resource_config(
        self, 
        partition_name: str,
        max_nodes: int,
        **overrides
    ) -> ResourceConfig:
        """
        Create a ResourceConfig from partition information.
        
        Args:
            partition_name: Name of the partition
            max_nodes: Maximum number of nodes
            **overrides: Additional ResourceConfig parameters to override
            
        Returns:
            ResourceConfig instance
        """
        partition = self.get_partition(partition_name)
        
        if not partition:
            logger.warning(f"Partition '{partition_name}' not found, using defaults")
            return ResourceConfig(max_nodes=max_nodes, **overrides)
        
        # Extract processor info
        processor = partition.processor or {}
        num_cpus = processor.get('num_cpus', 1)
        
        # Extract device info (GPUs)
        num_gpus = 0
        if partition.devices:
            num_gpus = partition.devices[0].get('num_devices', 0)
        
        # Parse access options for partition and account
        partition_opt = None
        account_opt = None
        for access in partition.access:
            if '--partition=' in access:
                partition_opt = access.split('=')[1]
            elif '--account=' in access:
                account_opt = access.split('=')[1]
        
        # Create ResourceConfig
        config = ResourceConfig(
            max_nodes=max_nodes,
            procs_per_node=num_cpus,
            gpus_per_node=num_gpus,
            partition=partition_opt,
            account=account_opt,
        )
        
        # Apply overrides
        for key, value in overrides.items():
            if hasattr(config, key):
                setattr(config, key, value)
        
        return config
    
    def create_backend_config(
        self,
        partition_name: str,
        **overrides
    ) -> BackendConfig:
        """
        Create a BackendConfig from partition information.
        
        Args:
            partition_name: Name of the partition
            **overrides: Additional BackendConfig parameters to override
            
        Returns:
            BackendConfig instance
        """
        partition = self.get_partition(partition_name)
        
        if not partition:
            logger.warning(f"Partition '{partition_name}' not found, using defaults")
            return BackendConfig(**overrides)
        
        # Map scheduler string to enum
        scheduler_map = {
            'local': SchedulerBackend.LOCAL,
            'slurm': SchedulerBackend.SLURM,
            'pbs': SchedulerBackend.PBS,
        }
        scheduler = scheduler_map.get(partition.scheduler, SchedulerBackend.LOCAL)
        
        # Map launcher string to enum (or keep as string for custom launchers)
        launcher_map = {
            'srun': LauncherBackend.SRUN,
            'mpirun': LauncherBackend.MPIRUN,
            'mpiexec': LauncherBackend.MPIEXEC,
            'simple': LauncherBackend.SIMPLE,
        }
        launcher = launcher_map.get(partition.launcher, partition.launcher)
        
        # Map module system string to enum
        module_map = {
            'nomod': ModuleBackend.NOMOD,
            'tmod': ModuleBackend.TMOD,
            'tmod4': ModuleBackend.TMOD4,
            'lmod': ModuleBackend.LMOD,
        }
        modules_system = ModuleBackend.NOMOD
        if self.current_system:
            modules_system = module_map.get(
                self.current_system.modules_system, 
                ModuleBackend.NOMOD
            )
        
        # Create BackendConfig
        config = BackendConfig(
            scheduler=scheduler,
            launcher=launcher,
            module_system=modules_system,
            scheduler_options=partition.sched_options,
        )
        
        # Apply overrides
        for key, value in overrides.items():
            if hasattr(config, key):
                setattr(config, key, value)
        
        return config
    
    def get_environment_modules(self, environment_name: str) -> list[str]:
        """
        Get the list of modules for a specific environment.
        
        Args:
            environment_name: Name of the environment
            
        Returns:
            List of module names
        """
        if not self.site_config:
            return []
        
        env = self.site_config.get_environment(environment_name)
        if env:
            return env.modules
        
        return []
    
    def print_system_info(self):
        """Print information about the current system configuration."""
        if not self.current_system:
            print("No system configuration loaded")
            return
        
        print(f"\nSystem: {self.current_system.name}")
        print(f"Description: {self.current_system.descr}")
        print(f"Module System: {self.current_system.modules_system}")
        print(f"\nPartitions:")
        
        for partition in self.current_system.partitions:
            print(f"\n  {partition.name}:")
            print(f"    Scheduler: {partition.scheduler}")
            print(f"    Launcher: {partition.launcher}")
            
            if partition.processor:
                print(f"    CPUs per node: {partition.processor.get('num_cpus', 'N/A')}")
            
            if partition.devices:
                gpu = partition.devices[0]
                print(f"    GPUs per node: {gpu.get('num_devices', 'N/A')}")
                print(f"    GPU Model: {gpu.get('model', 'N/A')}")


def auto_load_system_config(
    search_paths: Optional[list[Path]] = None
) -> Optional[SystemConfigLoader]:
    """
    Automatically search for and load a system configuration file.
    
    Args:
        search_paths: Optional list of paths to search. 
                     If None, searches current directory and parent.
    
    Returns:
        SystemConfigLoader if a configuration was found, None otherwise
    """
    if search_paths is None:
        search_paths = [
            Path.cwd(),
            Path.cwd().parent,
        ]
    
    loader = SystemConfigLoader()
    config_file = loader.find_config_file(search_paths)
    
    if config_file:
        if loader.load_configuration(config_file):
            return loader
    
    logger.info("No system configuration file found")
    return None
