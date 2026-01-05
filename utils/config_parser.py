"""YAML configuration parser for HPC-ScaleTest.

This module handles parsing of run.yaml configuration files and converts them
to orchestrator configuration objects.
"""

import logging
import yaml
from pathlib import Path
from typing import Optional, Dict, Any, List

from utils.validators import ConfigValidator, ValidationError
from utils.parsers import parse_tuple  # Use unified parser

logger = logging.getLogger(__name__)


class YAMLConfigParser:
    """Parser for YAML configuration files."""
    
    def __init__(self, config_file: Path):
        """
        Initialize YAML config parser.
        
        Args:
            config_file: Path to YAML configuration file
        """
        self.config_file = Path(config_file)
        self.config_data: Optional[Dict[str, Any]] = None
    
    def parse(self) -> Dict[str, Any]:
        """
        Parse YAML configuration file.
        
        Returns:
            Dictionary containing parsed configuration
            
        Raises:
            FileNotFoundError: If config file doesn't exist
            yaml.YAMLError: If config file has invalid YAML syntax
        """
        if not self.config_file.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.config_file}")
        
        logger.info(f"Parsing configuration from {self.config_file}")
        
        try:
            with open(self.config_file, 'r') as f:
                self.config_data = yaml.safe_load(f)
            
            if self.config_data is None:
                logger.warning("Configuration file is empty, using defaults")
                self.config_data = {}
            
            logger.info(f"Successfully parsed configuration")
            return self.config_data
            
        except yaml.YAMLError as e:
            logger.error(f"Failed to parse YAML configuration: {e}")
            raise
        except Exception as e:
            logger.error(f"Error reading configuration file: {e}")
            raise
    
    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value by key.
        
        Args:
            key: Configuration key (supports nested keys with dots, e.g., 'scaling.type')
            default: Default value if key not found
            
        Returns:
            Configuration value or default
        """
        if self.config_data is None:
            self.parse()
        
        # Handle nested keys
        keys = key.split('.')
        value = self.config_data
        
        for k in keys:
            if isinstance(value, dict) and k in value:
                value = value[k]
            else:
                return default
        
        return value
    
    def to_orchestrator_config(self) -> Dict[str, Any]:
        """
        Convert YAML configuration to orchestrator configuration format.
        
        Returns:
            Dictionary suitable for OrchestratorConfig initialization
        """
        if self.config_data is None:
            self.parse()
        
        config = {}
        
        # Source configuration (required)
        if 'repository' in self.config_data:
            config['source'] = self.config_data['repository']
        elif 'source' in self.config_data:
            config['source'] = self.config_data['source']
        else:
            raise ValueError("Configuration must specify 'repository' or 'source'")
        
        # Scaling configuration
        if 'scaling' in self.config_data:
            scaling = self.config_data['scaling']
            if isinstance(scaling, str):
                # Simple format: scaling: weak
                config['scaling_type'] = scaling
            else:
                # Detailed format with sub-keys
                config['scaling_type'] = scaling.get('type', 'strong')
                if 'nodes' in scaling:
                    config['max_nodes'] = scaling['nodes']
                if 'initial_procs' in scaling:
                    config['initial_procs'] = self._parse_tuple(scaling['initial_procs'])
                if 'initial_domain' in scaling:
                    config['initial_domain'] = self._parse_tuple(scaling['initial_domain'])
                if 'initial_cells' in scaling:
                    config['initial_cells'] = self._parse_tuple(scaling['initial_cells'])
                if 'particles_per_cell' in scaling:
                    config['particles_per_cell'] = self._parse_tuple(scaling['particles_per_cell'])
        
        # GENERIC: Variable mapping for input files
        if 'variable_map' in self.config_data:
            config['variable_map'] = self.config_data['variable_map']
        
        # Scaling factor (if defined, enables full weak scaling mode)
        if 'scaling_factor' in self.config_data:
            config['scaling_factor'] = float(self.config_data['scaling_factor'])
        
        # Scaling dimensions (2D or 3D cyclic pattern)
        if 'scaling_dimensions' in self.config_data:
            config['scaling_dimensions'] = int(self.config_data['scaling_dimensions'])
        
        # Scaling factors (weak_scaling_factors, strong_scaling_factors)
        if 'weak_scaling_factors' in self.config_data:
            factors = self.config_data['weak_scaling_factors']
            if isinstance(factors, list):
                config['weak_scaling_factors'] = [float(f) for f in factors]
        
        if 'strong_scaling_factors' in self.config_data:
            factors = self.config_data['strong_scaling_factors']
            if isinstance(factors, list):
                config['strong_scaling_factors'] = [float(f) for f in factors]
        
        # Node configuration (top-level) - supports both 'nodes' and 'max_nodes'
        if 'max_nodes' in self.config_data:
            config['max_nodes'] = int(self.config_data['max_nodes'])
        elif 'nodes' in self.config_data:
            config['max_nodes'] = int(self.config_data['nodes'])
        
        # Top-level initial_procs, initial_domain, initial_cells (alternative to nested in scaling)
        if 'initial_procs' in self.config_data:
            config['initial_procs'] = self._parse_tuple(self.config_data['initial_procs'])
        
        if 'initial_domain' in self.config_data:
            config['initial_domain'] = self._parse_tuple(self.config_data['initial_domain'])
        
        if 'initial_cells' in self.config_data:
            config['initial_cells'] = self._parse_tuple(self.config_data['initial_cells'])
        
        # Top-level particles_per_cell (alternative to nested in scaling)
        if 'particles_per_cell' in self.config_data:
            ppc = self.config_data['particles_per_cell']
            if isinstance(ppc, dict):
                # Format: particles_per_cell: {x: 2, y: 2, z: 2}
                config['particles_per_cell'] = (ppc.get('x', 1), ppc.get('y', 1), ppc.get('z', 1))
            else:
                # Format: particles_per_cell: [2, 2, 2]
                config['particles_per_cell'] = self._parse_tuple(ppc)
        
        # Hardware configuration
        # Support both 'hardware' nested format and 'hardware_type' flat format
        if 'hardware_type' in self.config_data:
            config['hardware_type'] = self.config_data['hardware_type']
        elif 'hardware' in self.config_data:
            hw = self.config_data['hardware']
            if isinstance(hw, str):
                config['hardware_type'] = hw
            else:
                config['hardware_type'] = hw.get('type', 'cpu')
                if 'procs_per_node' in hw:
                    config['procs_per_node'] = hw['procs_per_node']
                if 'cpus_per_node' in hw:
                    config['cpus_per_node'] = hw['cpus_per_node']
                if 'gpus_per_node' in hw:
                    config['gpus_per_node'] = hw['gpus_per_node']
        
        # GPU configuration (top-level) - used by GPU jobs
        if 'gpus_per_node' in self.config_data:
            config['gpus_per_node'] = int(self.config_data['gpus_per_node'])
        
        # Top-level procs_per_node (alternative to nested in hardware)
        if 'procs_per_node' in self.config_data:
            config['procs_per_node'] = int(self.config_data['procs_per_node'])
        
        # Top-level cpus_per_node (total CPU cores per node for SLURM allocation)
        if 'cpus_per_node' in self.config_data:
            config['cpus_per_node'] = int(self.config_data['cpus_per_node'])
        
        # Resource configuration
        if 'partition' in self.config_data:
            config['partition'] = self.config_data['partition']
        if 'qos' in self.config_data:
            config['qos'] = self.config_data['qos']
        if 'qos_mapping' in self.config_data:
            config['qos_mapping'] = self.config_data['qos_mapping']
        if 'account' in self.config_data:
            config['account'] = self.config_data['account']
        if 'time_limit' in self.config_data:
            config['time_limit'] = self.config_data['time_limit']
        
        # Backend configuration
        if 'scheduler' in self.config_data:
            config['scheduler'] = self.config_data['scheduler']
        if 'launcher' in self.config_data:
            config['launcher'] = self.config_data['launcher']
        if 'module_system' in self.config_data:
            config['module_system'] = self.config_data['module_system']
        
        # Build configuration
        if 'build_system' in self.config_data:
            config['build_system'] = self.config_data['build_system']
        if 'modules' in self.config_data:
            modules = self.config_data['modules']
            if isinstance(modules, list):
                config['modules'] = modules
            elif isinstance(modules, str):
                config['modules'] = [m.strip() for m in modules.split(',')]
        
        # Executable configuration
        if 'executable' in self.config_data:
            config['executable_name'] = self.config_data['executable']
        if 'args' in self.config_data:
            args = self.config_data['args']
            if isinstance(args, list):
                config['executable_args'] = args
            elif isinstance(args, str):
                config['executable_args'] = [args]
        
        # Build options (cmake flags, etc.)
        if 'build_options' in self.config_data:
            build_opts = self.config_data['build_options']
            if isinstance(build_opts, dict):
                if 'cmake_opts' in build_opts:
                    # Parse cmake options into a dict
                    cmake_opts = build_opts['cmake_opts']
                    if isinstance(cmake_opts, str):
                        # Parse "-DKEY=VALUE -DKEY2=VALUE2" format
                        flags = {}
                        for part in cmake_opts.split():
                            if part.startswith('-D') and '=' in part:
                                key_val = part[2:].split('=', 1)
                                if len(key_val) == 2:
                                    flags[key_val[0]] = key_val[1]
                        config['build_flags'] = flags
                    elif isinstance(cmake_opts, dict):
                        config['build_flags'] = cmake_opts
        
        # Input file configuration
        if 'input_file' in self.config_data:
            config['input_file'] = self.config_data['input_file']
        if 'input_name' in self.config_data:
            config['input_file_name'] = self.config_data['input_name']
        
        # Behavior configuration
        # Note: 'verbose' is handled by logging setup in main(), not OrchestratorConfig
        if 'auto_submit' in self.config_data:
            config['auto_submit_jobs'] = self.config_data['auto_submit']
        if 'cleanup' in self.config_data:
            config['cleanup_after_build'] = self.config_data['cleanup']
        if 'no_reports' in self.config_data:
            config['generate_reports'] = not self.config_data['no_reports']
        
        # Test name
        if 'test_name' in self.config_data:
            config['test_name'] = self.config_data['test_name']
        
        # Output directories
        if 'output_dir' in self.config_data:
            config['output_dir'] = Path(self.config_data['output_dir'])
        if 'workspace_dir' in self.config_data:
            config['workspace_dir'] = Path(self.config_data['workspace_dir'])
        
        # LLM configuration (for intelligent parameter mapping)
        if 'use_llm_discovery' in self.config_data:
            config['use_llm_discovery'] = self.config_data['use_llm_discovery']
        if 'openai_api_key' in self.config_data:
            config['openai_api_key'] = self.config_data['openai_api_key']
        if 'llm_model' in self.config_data:
            config['llm_model'] = self.config_data['llm_model']
        if 'parameter_mapping' in self.config_data:
            config['parameter_mapping'] = self.config_data['parameter_mapping']
        
        logger.debug(f"Converted YAML config to orchestrator config: {config}")
        return config
    
    def _parse_tuple(self, value: Any) -> tuple:
        """Parse a tuple value from YAML - delegates to unified parser."""
        from utils.parsers import parse_tuple
        return parse_tuple(value)
    
    def validate(self) -> bool:
        """
        Validate configuration file.
        
        Returns:
            True if configuration is valid
            
        Raises:
            ValueError: If required fields are missing or invalid
        """
        if self.config_data is None:
            self.parse()
        
        # Check required fields
        if 'repository' not in self.config_data and 'source' not in self.config_data:
            raise ValueError("Configuration must specify 'repository' or 'source'")
        
        # REQUIRED: max_nodes must be specified (no hardcoded defaults)
        has_max_nodes = (
            'max_nodes' in self.config_data or
            'nodes' in self.config_data or
            (isinstance(self.config_data.get('scaling'), dict) and 'nodes' in self.config_data['scaling'])
        )
        if not has_max_nodes:
            raise ValueError("Configuration MUST specify 'max_nodes' or 'nodes' - no hardcoded defaults allowed")
        
        # REQUIRED: partition must be specified for SLURM jobs
        if 'partition' not in self.config_data:
            raise ValueError("Configuration MUST specify 'partition' - no hardcoded defaults allowed")
        
        # Determine hardware type
        hw_type = 'cpu'
        if 'hardware_type' in self.config_data:
            hw_type = self.config_data['hardware_type']
        elif 'hardware' in self.config_data:
            hw = self.config_data['hardware']
            if isinstance(hw, str):
                hw_type = hw
            else:
                hw_type = hw.get('type', 'cpu')
        
        # GPU-specific required fields
        if hw_type.lower() == 'gpu':
            # gpus_per_node is REQUIRED for GPU jobs
            has_gpus_per_node = (
                'gpus_per_node' in self.config_data or
                (isinstance(self.config_data.get('hardware'), dict) and 'gpus_per_node' in self.config_data['hardware'])
            )
            if not has_gpus_per_node:
                raise ValueError("GPU configuration MUST specify 'gpus_per_node' - no hardcoded defaults allowed")
            
            # cpus_per_node or procs_per_node is REQUIRED for GPU jobs
            has_cpu_count = (
                'cpus_per_node' in self.config_data or
                'procs_per_node' in self.config_data or
                (isinstance(self.config_data.get('hardware'), dict) and 
                 ('cpus_per_node' in self.config_data['hardware'] or 'procs_per_node' in self.config_data['hardware']))
            )
            if not has_cpu_count:
                raise ValueError(
                    "GPU configuration MUST specify 'cpus_per_node' or 'procs_per_node' (total CPU cores per node, e.g., 32) - "
                    "this is required for proper SLURM allocation and MPI mapping"
                )
        
        # Validate scaling type if specified
        if 'scaling' in self.config_data:
            scaling = self.config_data['scaling']
            if isinstance(scaling, str):
                scaling_type = scaling
            else:
                scaling_type = scaling.get('type', 'strong')
            
            if scaling_type not in ['strong', 'weak']:
                raise ValueError(f"Invalid scaling type: {scaling_type}. Must be 'strong' or 'weak'")
        
        # Validate hardware type if specified
        if hw_type not in ['cpu', 'gpu']:
            raise ValueError(f"Invalid hardware type: {hw_type}. Must be 'cpu' or 'gpu'")
        
        logger.info("Configuration validation passed")
        return True


def load_yaml_config(config_file: Path) -> Dict[str, Any]:
    """
    Load and parse YAML configuration file with validation.
    
    Args:
        config_file: Path to YAML configuration file
        
    Returns:
        Dictionary containing parsed configuration for OrchestratorConfig
        
    Raises:
        ValidationError: If configuration is invalid
    """
    parser = YAMLConfigParser(config_file)
    parser.parse()
    
    # Validate configuration
    validator = ConfigValidator()
    try:
        validator.validate_and_raise(parser.config_data)
        logger.info("✓ Configuration validation passed")
    except ValidationError as e:
        logger.error(str(e))
        raise
    
    return parser.to_orchestrator_config()
