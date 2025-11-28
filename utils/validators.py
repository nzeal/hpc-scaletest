"""
Input validation for HPC-ScaleTest configurations.

Validates YAML configurations and simulation parameters to prevent runtime errors.
"""

import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple

logger = logging.getLogger(__name__)


class ValidationError(Exception):
    """Raised when validation fails."""
    pass


class ConfigValidator:
    """Validate YAML configuration files."""
    
    def validate(self, config: Dict[str, Any]) -> List[str]:
        """
        Validate configuration dictionary.
        
        Args:
            config: Configuration dictionary from YAML
            
        Returns:
            List of error messages (empty if valid)
        """
        errors = []
        
        # Required fields
        errors.extend(self._validate_required_fields(config))
        
        # Source path
        errors.extend(self._validate_source(config))
        
        # Scaling configuration
        errors.extend(self._validate_scaling(config))
        
        # Resource configuration
        errors.extend(self._validate_resources(config))
        
        # Input file
        errors.extend(self._validate_input_file(config))
        
        # LLM configuration
        errors.extend(self._validate_llm_config(config))
        
        return errors
    
    def _validate_required_fields(self, config: Dict[str, Any]) -> List[str]:
        """Validate required configuration fields."""
        errors = []
        
        if 'repository' not in config and 'source' not in config:
            errors.append("Missing required field: 'repository' or 'source'")
        
        return errors
    
    def _validate_source(self, config: Dict[str, Any]) -> List[str]:
        """Validate source path/URL."""
        errors = []
        
        source = config.get('repository') or config.get('source')
        if not source:
            return errors
        
        # Check if it's a URL or local path
        if not source.startswith(('http://', 'https://', 'git@')):
            # Local path - verify it exists
            path = Path(source)
            if not path.exists():
                errors.append(f"Source path does not exist: {source}")
            elif not path.is_dir():
                errors.append(f"Source must be a directory: {source}")
        
        return errors
    
    def _validate_scaling(self, config: Dict[str, Any]) -> List[str]:
        """Validate scaling configuration."""
        errors = []
        
        # Scaling type
        scaling_type = config.get('scaling', 'strong')
        if scaling_type not in ['strong', 'weak']:
            errors.append(f"Invalid scaling type: '{scaling_type}'. Must be 'strong' or 'weak'")
        
        # Node count
        nodes = config.get('nodes', config.get('max_nodes'))
        if nodes is not None:
            if not isinstance(nodes, int):
                errors.append(f"'nodes' must be an integer, got {type(nodes).__name__}")
            elif nodes < 1:
                errors.append(f"'nodes' must be >= 1, got {nodes}")
            elif nodes > 10000:
                errors.append(f"'nodes' seems unreasonably large: {nodes}. Please verify.")
        
        # Weak scaling specific
        if scaling_type == 'weak':
            if not config.get('initial_domain') and not config.get('initial_cells'):
                errors.append("Weak scaling requires 'initial_domain' or 'initial_cells'")
            
            # Validate domain format
            if 'initial_domain' in config:
                domain = config['initial_domain']
                if not isinstance(domain, (list, tuple)) or len(domain) != 3:
                    errors.append(f"'initial_domain' must be [Lx, Ly, Lz], got {domain}")
                elif not all(isinstance(x, (int, float)) and x > 0 for x in domain):
                    errors.append(f"'initial_domain' values must be positive numbers")
            
            # Validate cells format
            if 'initial_cells' in config:
                cells = config['initial_cells']
                if not isinstance(cells, (list, tuple)) or len(cells) != 3:
                    errors.append(f"'initial_cells' must be [nx, ny, nz], got {cells}")
                elif not all(isinstance(x, int) and x > 0 for x in cells):
                    errors.append(f"'initial_cells' values must be positive integers")
        
        return errors
    
    def _validate_resources(self, config: Dict[str, Any]) -> List[str]:
        """Validate resource configuration."""
        errors = []
        
        # Partition
        partition = config.get('partition')
        if partition == 'X_usr_prod':
            errors.append("Please specify actual partition name (default 'X_usr_prod' is placeholder)")
        
        # Account
        account = config.get('account')
        if account == 'cin_X':
            errors.append("Please specify actual account name (default 'cin_X' is placeholder)")
        
        # Processors per node
        procs = config.get('procs_per_node')
        if procs is not None:
            if not isinstance(procs, int):
                errors.append(f"'procs_per_node' must be integer, got {type(procs).__name__}")
            elif procs < 1:
                errors.append(f"'procs_per_node' must be >= 1, got {procs}")
            elif procs > 512:
                errors.append(f"'procs_per_node' seems unreasonably large: {procs}")
        
        # GPUs per node
        gpus = config.get('gpus_per_node')
        if gpus is not None and gpus < 0:
            errors.append(f"'gpus_per_node' cannot be negative, got {gpus}")
        
        # Hardware type
        hardware = config.get('hardware', config.get('hardware_type', 'cpu'))
        if hardware not in ['cpu', 'gpu']:
            errors.append(f"Invalid hardware type: '{hardware}'. Must be 'cpu' or 'gpu'")
        
        return errors
    
    def _validate_input_file(self, config: Dict[str, Any]) -> List[str]:
        """Validate input file configuration."""
        errors = []
        
        input_file = config.get('input_file')
        if input_file:
            path = Path(input_file)
            # Only validate if it's not a relative path inside the source
            if path.is_absolute() and not path.exists():
                errors.append(f"Input file not found: {input_file}")
        
        return errors
    
    def _validate_llm_config(self, config: Dict[str, Any]) -> List[str]:
        """Validate LLM configuration."""
        errors = []
        
        use_llm = config.get('use_llm_discovery', False)
        if use_llm:
            # Check for API key
            api_key = config.get('openai_api_key')
            if not api_key:
                import os
                if not os.getenv('OPENAI_API_KEY'):
                    errors.append(
                        "LLM discovery enabled but no API key found. "
                        "Set 'openai_api_key' in YAML or OPENAI_API_KEY env var"
                    )
        
        # Validate parameter mapping format
        param_map = config.get('parameter_mapping')
        if param_map:
            if not isinstance(param_map, dict):
                errors.append(f"'parameter_mapping' must be a dictionary, got {type(param_map).__name__}")
            else:
                for key, value in param_map.items():
                    if not isinstance(value, list):
                        errors.append(f"parameter_mapping['{key}'] must be a list, got {type(value).__name__}")
        
        return errors
    
    def validate_and_raise(self, config: Dict[str, Any]):
        """
        Validate configuration and raise ValidationError if invalid.
        
        Args:
            config: Configuration dictionary
            
        Raises:
            ValidationError: If validation fails
        """
        errors = self.validate(config)
        if errors:
            error_msg = "\n❌ Configuration Validation Failed:\n" + "\n".join(f"  • {err}" for err in errors)
            raise ValidationError(error_msg)


class ParameterValidator:
    """Validate simulation parameters."""
    
    def validate_decomposition(
        self,
        nx: int,
        ny: int,
        nz: int,
        decomp_x: int,
        decomp_y: int,
        decomp_z: int
    ) -> List[str]:
        """
        Validate domain decomposition against grid dimensions.
        
        This prevents the critical "ZLEN=28 doesn't divide nz=1" error.
        
        Args:
            nx, ny, nz: Grid dimensions
            decomp_x, decomp_y, decomp_z: MPI decomposition
            
        Returns:
            List of error messages (empty if valid)
        """
        errors = []
        
        # Check divisibility
        if nx % decomp_x != 0:
            errors.append(
                f"Grid X ({nx}) not divisible by DecompX ({decomp_x}). "
                f"Remainder: {nx % decomp_x}"
            )
        
        if ny % decomp_y != 0:
            errors.append(
                f"Grid Y ({ny}) not divisible by DecompY ({decomp_y}). "
                f"Remainder: {ny % decomp_y}"
            )
        
        if nz % decomp_z != 0:
            errors.append(
                f"Grid Z ({nz}) not divisible by DecompZ ({decomp_z}). "
                f"Remainder: {nz % decomp_z}"
            )
        
        # Critical check: Never split single-cell dimensions
        if nz == 1 and decomp_z > 1:
            errors.append(
                f"CRITICAL: Cannot split nz=1 into {decomp_z} subdomains. "
                f"This will cause 'bad_array_new_length' error. Set DecompZ=1."
            )
        
        if nx == 1 and decomp_x > 1:
            errors.append(f"Cannot split nx=1 into {decomp_x} subdomains")
        
        if ny == 1 and decomp_y > 1:
            errors.append(f"Cannot split ny=1 into {decomp_y} subdomains")
        
        # Check for very small subdomains (warning)
        cells_x = nx / decomp_x if decomp_x > 0 else 0
        cells_y = ny / decomp_y if decomp_y > 0 else 0
        cells_z = nz / decomp_z if decomp_z > 0 else 0
        
        if cells_x < 2:
            errors.append(
                f"WARNING: Very small X subdomain ({cells_x:.1f} cells). "
                f"Consider reducing DecompX or increasing nx"
            )
        
        if cells_y < 2:
            errors.append(
                f"WARNING: Very small Y subdomain ({cells_y:.1f} cells). "
                f"Consider reducing DecompY or increasing ny"
            )
        
        if cells_z < 2 and nz > 1:
            errors.append(
                f"WARNING: Very small Z subdomain ({cells_z:.1f} cells). "
                f"Consider reducing DecompZ or increasing nz"
            )
        
        return errors
    
    def validate_and_raise(
        self,
        nx: int,
        ny: int,
        nz: int,
        decomp_x: int,
        decomp_y: int,
        decomp_z: int
    ):
        """
        Validate and raise ValidationError if invalid.
        
        Raises:
            ValidationError: If decomposition is invalid
        """
        errors = self.validate_decomposition(nx, ny, nz, decomp_x, decomp_y, decomp_z)
        if errors:
            error_msg = "\n❌ Decomposition Validation Failed:\n" + "\n".join(f"  • {err}" for err in errors)
            raise ValidationError(error_msg)


def validate_job_config(job_config, test) -> Tuple[bool, List[str]]:
    """
    Pre-execution validation of job configuration.
    
    Prevents runtime crashes by catching configuration errors before submission.
    
    Args:
        job_config: JobConfig object with scaled grid dimensions
        test: Test object with simulation parameters
        
    Returns:
        Tuple of (is_valid, error_messages)
    """
    errors = []
    
    # Validate decomposition using SCALED grid dimensions from job_config, NOT initial_cells
    # This is CRITICAL: Each job has different grid dimensions (weak scaling)
    if job_config.cell_count is not None:
        nx, ny, nz = job_config.cell_count
        px, py, pz = job_config.procs_decomposition
        
        validator = ParameterValidator()
        decomp_errors = validator.validate_decomposition(nx, ny, nz, px, py, pz)
        errors.extend(decomp_errors)
    elif hasattr(test.scaling_config, 'initial_cells') and test.scaling_config.initial_cells:
        # Fallback to initial_cells if cell_count not set (should not happen)
        nx, ny, nz = test.scaling_config.initial_cells
        px, py, pz = job_config.procs_decomposition
        
        validator = ParameterValidator()
        decomp_errors = validator.validate_decomposition(nx, ny, nz, px, py, pz)
        errors.extend(decomp_errors)
    
    # Validate total process count
    expected_procs = (job_config.procs_decomposition[0] *
                     job_config.procs_decomposition[1] *
                     job_config.procs_decomposition[2])
    
    if expected_procs != job_config.num_procs:
        errors.append(
            f"Process count mismatch: DecompX×DecompY×DecompZ = {expected_procs}, "
            f"but num_procs = {job_config.num_procs}"
        )
    
    is_valid = len(errors) == 0
    
    if not is_valid:
        logger.error("Job configuration validation failed:")
        for error in errors:
            logger.error(f"  • {error}")
    
    return is_valid, errors
