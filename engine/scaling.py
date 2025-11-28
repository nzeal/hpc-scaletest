#!/usr/bin/env python3
"""
scaling.py - 1D/2D/3D weak scaling engine for PIC simulation inputs
with strict rules for input file parsing and deterministic weak scaling.

This implementation follows strict rules for:
- Node 1: Always override with values from run.yaml, ignore input file values
- Node > 1: Apply strict weak-scaling rules that preserve previous node values
- MPI topology correction for Node 4+ becomes correct (preserving previous values)
"""

import os
import yaml
import copy
import argparse
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# ------------------- Colors -------------------
class Colors:
    """ANSI color codes for terminal output."""
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    RESET = '\033[0m'

# ------------------- Safe conversions -------------------
def safe_float(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return None

def safe_int(value):
    try:
        return int(value)
    except (ValueError, TypeError):
        return None

def scale_value(value, factor):
    num = safe_float(value)
    if num is not None:
        return num * factor
    # If value is not a valid float, return 1.0 (safe default for scaling)
    return 1.0

def scale_int(value, factor):
    num = safe_int(value)
    if num is not None:
        return int(round(num * factor))
    # If value is not a valid integer, return 1 (safe default for scaling)
    return 1

# ------------------- Input File Parser -------------------
class InputFileParser:
    """Parse and extract parameters from PIC input files."""
    
    def __init__(self):
        # Map of parameter names to extract
        self.param_names = [
            'Lx', 'Ly', 'Lz',
            'nxc', 'nyc', 'nzc',
            'XLEN', 'YLEN', 'ZLEN',
            'npcelx', 'npcely', 'npcelz'
        ]
    
    def parse_input_file(self, input_file):
        """Parse input file and extract parameter values."""
        if not os.path.exists(input_file):
            logger.warning(f"Input file not found: {input_file}")
            return {}
        
        with open(input_file, 'r') as f:
            content = f.read()
        
        params = {}
        for param in self.param_names:
            value = self._extract_parameter(content, param)
            if value is not None:
                params[param] = value
        
        return params
    
    def _extract_parameter(self, content, param_name):
        """Extract parameter value from content."""
        # Look for patterns like "param = value" or "param: value"
        import re
        patterns = [
            rf'{re.escape(param_name)}\s*[=:]\s*([^\s#;\n]+)',
            rf'{re.escape(param_name)}\s+([^\s#;\n]+)'
        ]
        
        for pattern in patterns:
            match = re.search(pattern, content, re.IGNORECASE | re.MULTILINE)
            if match:
                value_str = match.group(1).strip()
                # Try to convert to numeric
                return self._safe_convert_to_numeric(value_str)
        
        return None
    
    def _safe_convert_to_numeric(self, value_str):
        """Safely convert string to numeric type."""
        if not value_str or value_str.isspace():
            return value_str
        
        try:
            if '.' in value_str or 'e' in value_str.lower():
                return float(value_str)
            else:
                return int(value_str)
        except ValueError:
            # Not numeric - return as string (could be placeholder)
            return value_str
    
    def generate_scaled_input(self, base_input, output_file, config):
        """Generate scaled input file with new parameter values."""
        if not os.path.exists(base_input):
            logger.error(f"Base input file not found: {base_input}")
            return False
        
        with open(base_input, 'r') as f:
            content = f.read()
        
        # Replace parameters with new values
        replacements = {
            'Lx': config.get('Lx'),
            'Ly': config.get('Ly'),
            'Lz': config.get('Lz'),
            'nxc': config.get('nx'),
            'nyc': config.get('ny'),
            'nzc': config.get('nz'),
            'XLEN': config.get('nproc_x'),
            'YLEN': config.get('nproc_y'),
            'ZLEN': config.get('nproc_z'),
            'npcelx': config.get('num_particles_x'),
            'npcely': config.get('num_particles_y'),
            'npcelz': config.get('num_particles_z')
        }
        
        # Also handle common typos in the input file
        typo_fixes = {
            'nyz': config.get('nz')  # Fix for nzc = nyz typo
        }
        
        # Apply main replacements
        for param, new_value in replacements.items():
            if new_value is not None:
                content = self._replace_parameter(content, param, new_value)
        
        # Apply typo fixes
        for param, new_value in typo_fixes.items():
            if new_value is not None:
                content = self._replace_parameter(content, param, new_value)
        
        # Write output file
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, 'w') as f:
            f.write(content)
        
        logger.info(f"Generated scaled input file: {output_file}")
        return True
    
    def _replace_parameter(self, content, param_name, new_value):
        """Replace parameter value in content."""
        import re
        # Pattern to match "param = value" or "param: value" with flexible spacing
        pattern = rf'({re.escape(param_name)}\s*[=:])\s*[^\s#;\n]+'
        replacement = f"\\g<1> {new_value}"
        
        # Try first pattern
        new_content = re.sub(pattern, replacement, content, flags=re.IGNORECASE | re.MULTILINE)
        if new_content != content:
            return new_content
        
        # Pattern to match "param value" (space-separated)
        pattern = rf'({re.escape(param_name)}\s+)\S+'
        replacement = f"\\g<1>{new_value}"
        new_content = re.sub(pattern, replacement, content, flags=re.IGNORECASE | re.MULTILINE)
        
        return new_content

# ------------------- Scaling Engine -------------------
class ScalingEngine:
    def __init__(self, input_file, output_dir, nodes, procs_per_node, scale_factor=2, dims=3, 
                 initial_domain=None, initial_cells=None, initial_procs=None, 
                 particles_per_cell=None, scaling_type="weak"):
        logger.info(f"Initializing Scaling Engine")
        logger.info(f"  Input file: {input_file}")
        logger.info(f"  Output directory: {output_dir}")
        logger.info(f"  Node sequence: {nodes}")
        logger.info(f"  Procs per node: {procs_per_node}")
        logger.info(f"  Scale factor: {scale_factor}")
        logger.info(f"  Scaling dimensions: {dims}D")
        logger.info(f"  Scaling type: {scaling_type}")
        
        self.input_file = input_file
        self.output_dir = output_dir
        self.nodes = nodes
        self.procs_per_node = procs_per_node
        self.scale_factor = scale_factor
        self.dims = dims  # 1, 2, or 3
        self.initial_domain = initial_domain
        self.initial_cells = initial_cells
        self.initial_procs = initial_procs
        self.particles_per_cell = particles_per_cell
        self.scaling_type = scaling_type.lower()
        
        # Flag to indicate if this is strong scaling
        self._is_strong_scaling = self.scaling_type == "strong"
        
        # Parse input file to get baseline parameters
        self.input_parser = InputFileParser()
        self.input_params = self.input_parser.parse_input_file(input_file)
        logger.info(f"  Parsed input parameters: {list(self.input_params.keys())}")
        
        # Create base configuration from YAML values (Node 1 override)
        # Initialize with sensible defaults to ensure keys always exist for scaling operations
        self.base_config = {
            'Lx': None,
            'Ly': None,
            'Lz': None,
            'nx': 1,      # Default cell count (will scale up)
            'ny': 1,      # Default cell count (will scale up)
            'nz': 1,      # Default cell count (will scale up)
            'nproc_x': 1, # Default MPI decomposition
            'nproc_y': 1, # Default MPI decomposition
            'nproc_z': 1, # Default MPI decomposition
            'num_particles_x': 1,
            'num_particles_y': 1,
            'num_particles_z': 1
        }
        
        # Override with initial values if provided (Node 1 uses these exactly)
        if initial_domain:
            self.base_config['Lx'] = initial_domain[0]
            self.base_config['Ly'] = initial_domain[1]
            self.base_config['Lz'] = initial_domain[2]
            logger.info(f"  Overriding domain with initial values: {initial_domain}")
        
        if initial_cells:
            self.base_config['nx'] = initial_cells[0]  # nx replaces nxc
            self.base_config['ny'] = initial_cells[1]  # ny replaces nyc
            self.base_config['nz'] = initial_cells[2]  # nz replaces nzc
            logger.info(f"  Overriding cells with initial values: {initial_cells}")
        
        if initial_procs:
            self.base_config['nproc_x'] = initial_procs[0]  # nproc_x replaces XLEN
            self.base_config['nproc_y'] = initial_procs[1]  # nproc_y replaces YLEN
            self.base_config['nproc_z'] = initial_procs[2]  # nproc_z replaces ZLEN
            logger.info(f"  Overriding MPI decomposition with initial values: {initial_procs}")
        
        # Set particles per cell (constant across all nodes)
        if particles_per_cell:
            self.base_config['num_particles_x'] = particles_per_cell[0]
            self.base_config['num_particles_y'] = particles_per_cell[1]
            self.base_config['num_particles_z'] = particles_per_cell[2]
            logger.info(f"  Setting particles per cell: {particles_per_cell}")
        
        # If any values are missing, try to get them from parsed input file
        # But only for non-critical parameters (Node 1 must use YAML values for critical ones)
        if self.base_config['Lx'] is None and 'Lx' in self.input_params:
            self.base_config['Lx'] = self.input_params['Lx']
        if self.base_config['Ly'] is None and 'Ly' in self.input_params:
            self.base_config['Ly'] = self.input_params['Ly']
        if self.base_config['Lz'] is None and 'Lz' in self.input_params:
            self.base_config['Lz'] = self.input_params['Lz']
        
        # Also try to get cell counts and decomposition from parsed input if not overridden
        if self.base_config['nx'] == 1 and 'nxc' in self.input_params:
            self.base_config['nx'] = self.input_params['nxc']
        if self.base_config['ny'] == 1 and 'nyc' in self.input_params:
            self.base_config['ny'] = self.input_params['nyc']
        if self.base_config['nz'] == 1 and 'nzc' in self.input_params:
            self.base_config['nz'] = self.input_params['nzc']
        
        # Get MPI decomposition from input if not overridden
        if self.base_config['nproc_x'] == 1 and 'XLEN' in self.input_params:
            self.base_config['nproc_x'] = self.input_params['XLEN']
        if self.base_config['nproc_y'] == 1 and 'YLEN' in self.input_params:
            self.base_config['nproc_y'] = self.input_params['YLEN']
        if self.base_config['nproc_z'] == 1 and 'ZLEN' in self.input_params:
            self.base_config['nproc_z'] = self.input_params['ZLEN']
        
        # Store all node configs for table
        self.all_configs = []
        
        # Output basename
        self.output_basename = os.path.splitext(os.path.basename(input_file))[0]
        os.makedirs(output_dir, exist_ok=True)

    def generate_all_cases(self):
        """Generate all scaling cases and write input files."""
        logger.info("\nGenerating scaling configurations for all nodes...")
        logger.info("=" * 80)
        
        for idx, node_count in enumerate(self.nodes):
            logger.info(f"\nProcessing Step {idx}: Node {node_count}")
            
            cfg = self.scale_node(idx, node_count)
            self.all_configs.append((node_count, cfg))
            
            # Generate input file for this node
            output_file = os.path.join(self.output_dir, f"{self.output_basename}_node{node_count}.in")
            success = self.input_parser.generate_scaled_input(self.input_file, output_file, cfg)
            if success:
                logger.info(f"  ✓ Generated: {output_file}")
            else:
                logger.error(f"  ✗ Failed to generate: {output_file}")
        
        logger.info("\n" + "=" * 80)
        logger.info("All configurations generated successfully\n")
        
        # Print color-coded verification table
        self.print_scaling_table()

    def scale_node(self, step_index, node_count):
        """Scale configuration for a specific node."""
        if step_index == 0:
            # Node 1: Baseline - use exact values from YAML, no scaling
            cfg = copy.deepcopy(self.base_config)
            logger.info(f"  {Colors.GREEN}BASELINE{Colors.RESET} - Using exact values from YAML")
            logger.info(f"    Domain: Lx={cfg.get('Lx', 'N/A')}, Ly={cfg.get('Ly', 'N/A')}, Lz={cfg.get('Lz', 'N/A')}")
            logger.info(f"    Cells:  nx={cfg.get('nx', 'N/A')}, ny={cfg.get('ny', 'N/A')}, nz={cfg.get('nz', 'N/A')}")
            logger.info(f"    MPI:    px={cfg.get('nproc_x', 'N/A')}, py={cfg.get('nproc_y', 'N/A')}, pz={cfg.get('nproc_z', 'N/A')}")
            cfg['active_dim'] = 'BASELINE'
            # Sanitize config to ensure all numeric values are proper types (not strings from input file)
            cfg = self._sanitize_config(cfg)
            return cfg
        
        # Check if this is strong scaling
        if self._is_strong_scaling:
            # For strong scaling, keep problem size constant and only adjust MPI decomposition
            cfg = copy.deepcopy(self.base_config)
            
            # Calculate total required processes
            required_ranks = node_count * self.procs_per_node
            
            # For strong scaling, we need to redistribute the same problem size across more processors
            # We'll adjust the MPI decomposition while keeping domain and cell counts constant
            logger.info(f"  {Colors.RED}Strong Scaling{Colors.RESET} - Keeping problem size constant, adjusting MPI decomposition")
            logger.info(f"    Required ranks: {required_ranks}")
            logger.info(f"    Domain: Lx={cfg.get('Lx', 'N/A')}, Ly={cfg.get('Ly', 'N/A')}, Lz={cfg.get('Lz', 'N/A')} (CONSTANT)")
            logger.info(f"    Cells:  nx={cfg.get('nx', 'N/A')}, ny={cfg.get('ny', 'N/A')}, nz={cfg.get('nz', 'N/A')} (CONSTANT)")
            
            # Adjust MPI decomposition to match required ranks while maintaining problem size
            cfg = self.correct_mpi_for_strong_scaling(cfg, required_ranks, node_count)
            cfg['active_dim'] = 'STRONG_SCALING'
            
            # Ensure all numeric fields are properly converted (not strings from input file)
            cfg = self._sanitize_config(cfg)
            
            return cfg
        else:
            # Node > 1: Apply weak scaling rules
            # Start with previous node's configuration (not base config)
            prev_node_count, prev_cfg = self.all_configs[-1]
            cfg = copy.deepcopy(prev_cfg)
            
            # Determine active dimension
            active_dim = self.get_active_dimension(step_index)
            factor = self.scale_factor
            
            logger.info(f"  {Colors.RED}Scaling dimension: {active_dim}{Colors.RESET} (factor={factor})")
            
            # Apply scaling to active dimension only
            if active_dim == 'X':
                cfg['Lx'] = scale_value(cfg.get('Lx'), factor)
                cfg['nx'] = scale_int(cfg.get('nx', 1), factor)
                cfg['nproc_x'] = scale_int(cfg.get('nproc_x', 1), factor)
            elif active_dim == 'Y':
                cfg['Ly'] = scale_value(cfg.get('Ly'), factor)
                cfg['ny'] = scale_int(cfg.get('ny', 1), factor)
                cfg['nproc_y'] = scale_int(cfg.get('nproc_y', 1), factor)
            elif active_dim == 'Z':
                cfg['Lz'] = scale_value(cfg.get('Lz'), factor)
                cfg['nz'] = scale_int(cfg.get('nz', 1), factor)
                cfg['nproc_z'] = scale_int(cfg.get('nproc_z', 1), factor)
            
            # Log scaling changes
            logger.info(f"    Domain: Lx={prev_cfg.get('Lx', 'N/A')}→{cfg.get('Lx', 'N/A')}, Ly={prev_cfg.get('Ly', 'N/A')}→{cfg.get('Ly', 'N/A')}, Lz={prev_cfg.get('Lz', 'N/A')}→{cfg.get('Lz', 'N/A')}")
            logger.info(f"    Cells:  nx={prev_cfg.get('nx', 'N/A')}→{cfg.get('nx', 'N/A')}, ny={prev_cfg.get('ny', 'N/A')}→{cfg.get('ny', 'N/A')}, nz={prev_cfg.get('nz', 'N/A')}→{cfg.get('nz', 'N/A')}")
            logger.info(f"    MPI:    px={prev_cfg.get('nproc_x', 'N/A')}→{cfg.get('nproc_x', 'N/A')}, py={prev_cfg.get('nproc_y', 'N/A')}→{cfg.get('nproc_y', 'N/A')}, pz={prev_cfg.get('nproc_z', 'N/A')}→{cfg.get('nproc_z', 'N/A')}")
            
            # MPI decomposition correction
            cfg = self.correct_mpi(cfg, active_dim, node_count)
            cfg['active_dim'] = active_dim
            
            # Ensure all numeric fields are properly converted (not strings from input file)
            cfg = self._sanitize_config(cfg)
            
            return cfg
    
    def _sanitize_config(self, cfg):
        """Convert string values to proper numeric types."""
        # Numeric integer fields
        int_fields = ['nx', 'ny', 'nz', 'nproc_x', 'nproc_y', 'nproc_z', 
                      'num_particles_x', 'num_particles_y', 'num_particles_z']
        for field in int_fields:
            if field in cfg:
                val = safe_int(cfg[field])
                cfg[field] = val if val is not None else 1
        
        # Numeric float fields
        float_fields = ['Lx', 'Ly', 'Lz']
        for field in float_fields:
            if field in cfg:
                val = safe_float(cfg[field])
                cfg[field] = val if val is not None else 1.0
        
        return cfg

    def get_active_dimension(self, step_index):
        """Determine active dimension based on scaling dimensions and step index."""
        # 1D always X
        if self.dims == 1:
            return 'X'
        # 2D alternating X/Y
        elif self.dims == 2:
            return 'X' if step_index % 2 == 1 else 'Y'
        # 3D cycling X/Y/Z
        elif self.dims == 3:
            cycle = ['X', 'Y', 'Z']
            return cycle[(step_index - 1) % 3]

    def correct_mpi(self, cfg, active_dim, node_count):
        """
        Correct MPI decomposition to match required ranks while ensuring divisibility with grid.
        
        CRITICAL: MPI decomposition must divide evenly into grid dimensions.
        When a grid dimension is not divisible by required MPI count, we SCALE THE GRID
        to the nearest size that IS divisible, maintaining weak scaling invariant.
        """
        required_ranks = node_count * self.procs_per_node
        
        # Handle None values and non-numeric values for MPI process counts
        nproc_x = cfg.get('nproc_x')
        nproc_y = cfg.get('nproc_y')
        nproc_z = cfg.get('nproc_z')
        
        # Convert to int, handling strings that aren't numbers
        nproc_x = safe_int(nproc_x) or 1
        nproc_y = safe_int(nproc_y) or 1
        nproc_z = safe_int(nproc_z) or 1
        
        # Update cfg with sanitized values
        cfg['nproc_x'] = nproc_x
        cfg['nproc_y'] = nproc_y
        cfg['nproc_z'] = nproc_z
        
        # Get grid dimensions (cells)
        nx = safe_int(cfg.get('nx')) or 1
        ny = safe_int(cfg.get('ny')) or 1
        nz = safe_int(cfg.get('nz')) or 1
        
        # Get domain size (physical)
        Lx = safe_float(cfg.get('Lx')) or 1.0
        Ly = safe_float(cfg.get('Ly')) or 1.0
        Lz = safe_float(cfg.get('Lz')) or 1.0
        
        actual_ranks = nproc_x * nproc_y * nproc_z
        
        if actual_ranks == required_ranks:
            # Rank count matches. But verify divisibility; if needed, scale grid.
            self._ensure_divisibility_with_grid_scaling(cfg, required_ranks)
            return cfg
        
        logger.info(f"  {Colors.YELLOW}MPI correction needed{Colors.RESET}: {actual_ranks} → {required_ranks} ranks")
        
        # Correct active dimension, ensuring divisibility by scaling grid if necessary
        if active_dim == 'X':
            nproc_x_new = required_ranks // (nproc_y * nproc_z)
            
            # Ensure nproc_x_new divides into grid nx (cells)
            # If not, scale nx to the nearest multiple of nproc_x_new
            if nx % nproc_x_new != 0:
                nx_old = nx
                nx_new = self._scale_to_multiple(nx, nproc_x_new)
                logger.warning(f"    Grid X not divisible by nproc_x={nproc_x_new}. Scaling: {nx_old} → {nx_new}")
                # Also scale Lx proportionally
                Lx_old = Lx
                Lx_new = Lx * (nx_new / nx_old) if nx_old > 0 else Lx
                cfg['nx'] = nx_new
                cfg['Lx'] = Lx_new
                logger.info(f"    Domain X scaled: {Lx_old} → {Lx_new}")
            
            if nproc_x_new * nproc_y * nproc_z != required_ranks:
                raise ValueError(f"MPI correction failed for node {node_count}: Cannot distribute {required_ranks} ranks")
            logger.info(f"    Adjusting nproc_x: {nproc_x} → {nproc_x_new} (grid_x={cfg.get('nx')}, divisible: {cfg.get('nx') % nproc_x_new == 0})")
            cfg['nproc_x'] = nproc_x_new
            
        elif active_dim == 'Y':
            nproc_y_new = required_ranks // (nproc_x * nproc_z)
            
            # Ensure nproc_y_new divides into grid ny
            if ny % nproc_y_new != 0:
                ny_old = ny
                ny_new = self._scale_to_multiple(ny, nproc_y_new)
                logger.warning(f"    Grid Y not divisible by nproc_y={nproc_y_new}. Scaling: {ny_old} → {ny_new}")
                # Also scale Ly proportionally
                Ly_old = Ly
                Ly_new = Ly * (ny_new / ny_old) if ny_old > 0 else Ly
                cfg['ny'] = ny_new
                cfg['Ly'] = Ly_new
                logger.info(f"    Domain Y scaled: {Ly_old} → {Ly_new}")
            
            if nproc_x * nproc_y_new * nproc_z != required_ranks:
                raise ValueError(f"MPI correction failed for node {node_count}: Cannot distribute {required_ranks} ranks")
            logger.info(f"    Adjusting nproc_y: {nproc_y} → {nproc_y_new} (grid_y={cfg.get('ny')}, divisible: {cfg.get('ny') % nproc_y_new == 0})")
            cfg['nproc_y'] = nproc_y_new
            
        elif active_dim == 'Z':
            nproc_z_new = required_ranks // (nproc_x * nproc_y)
            
            # Ensure nproc_z_new divides into grid nz
            if nz % nproc_z_new != 0:
                nz_old = nz
                nz_new = self._scale_to_multiple(nz, nproc_z_new)
                logger.warning(f"    Grid Z not divisible by nproc_z={nproc_z_new}. Scaling: {nz_old} → {nz_new}")
                # Also scale Lz proportionally
                Lz_old = Lz
                Lz_new = Lz * (nz_new / nz_old) if nz_old > 0 else Lz
                cfg['nz'] = nz_new
                cfg['Lz'] = Lz_new
                logger.info(f"    Domain Z scaled: {Lz_old} → {Lz_new}")
            
            if nproc_x * nproc_y * nproc_z_new != required_ranks:
                raise ValueError(f"MPI correction failed for node {node_count}: Cannot distribute {required_ranks} ranks")
            logger.info(f"    Adjusting nproc_z: {nproc_z} → {nproc_z_new} (grid_z={cfg.get('nz')}, divisible: {cfg.get('nz') % nproc_z_new == 0})")
            cfg['nproc_z'] = nproc_z_new
        
        return cfg
    
    def correct_mpi_for_strong_scaling(self, cfg, required_ranks, node_count):
        """
        Correct MPI decomposition for strong scaling to match required ranks
        while keeping the problem size constant.
        
        For strong scaling, we redistribute the same problem across more processors
        by adjusting the MPI decomposition. CRITICAL: MPI decomposition must divide
        evenly into grid dimensions to avoid "bad_array_new_length" errors.
        
        Args:
            cfg: Configuration dictionary
            required_ranks: Total number of MPI ranks needed
            node_count: Number of nodes being used
            
        Returns:
            Updated configuration with corrected MPI decomposition
        """
        # Get current grid dimensions (fixed for strong scaling)
        nx = safe_int(cfg.get('nx')) or 1
        ny = safe_int(cfg.get('ny')) or 1
        nz = safe_int(cfg.get('nz')) or 1
        
        # Get current MPI decomposition
        nproc_x = safe_int(cfg.get('nproc_x')) or 1
        nproc_y = safe_int(cfg.get('nproc_y')) or 1
        nproc_z = safe_int(cfg.get('nproc_z')) or 1
        
        logger.info(f"    Current MPI decomposition: {nproc_x} × {nproc_y} × {nproc_z} = {nproc_x * nproc_y * nproc_z}")
        logger.info(f"    Required ranks: {required_ranks}")
        logger.info(f"    Grid dimensions: nx={nx}, ny={ny}, nz={nz}")
        
        # For strong scaling with fixed grid, we need MPI decomposition to divide evenly into grid
        # Get all divisors of required_ranks and grid dimensions
        new_x, new_y, new_z = self._find_compatible_decomposition(required_ranks, nx, ny, nz)
        
        logger.info(f"    Compatible MPI decomposition: {new_x} × {new_y} × {new_z} = {new_x * new_y * new_z}")
        
        # Verify divisibility
        if nx % new_x != 0 or ny % new_y != 0 or nz % new_z != 0:
            logger.warning(f"    Warning: Found decomposition may not be perfectly divisible with grid")
            logger.warning(f"      nx % new_x = {nx} % {new_x} = {nx % new_x if new_x > 0 else 'N/A'}")
            logger.warning(f"      ny % new_y = {ny} % {new_y} = {ny % new_y if new_y > 0 else 'N/A'}")
            logger.warning(f"      nz % new_z = {nz} % {new_z} = {nz % new_z if new_z > 0 else 'N/A'}")
        
        # Update configuration
        cfg['nproc_x'] = new_x
        cfg['nproc_y'] = new_y
        cfg['nproc_z'] = new_z
        
        return cfg
    
    def _find_compatible_decomposition(self, required_ranks, nx, ny, nz):
        """
        Find MPI decomposition that:
        1. Multiplies to required_ranks
        2. Divides evenly into grid dimensions (nx, ny, nz)
        
        This is critical for strong scaling to avoid "bad_array_new_length" errors.
        
        Args:
            required_ranks: Total MPI ranks needed
            nx, ny, nz: Grid dimensions (fixed for strong scaling)
            
        Returns:
            Tuple (nproc_x, nproc_y, nproc_z) that satisfies both constraints
        """
        # Get divisors of required_ranks
        rank_divisors = self._get_divisors(required_ranks)
        logger.debug(f"    Divisors of required_ranks {required_ranks}: {rank_divisors}")
        
        # Get divisors of grid dimensions
        x_divisors = self._get_divisors(nx)
        y_divisors = self._get_divisors(ny)
        z_divisors = self._get_divisors(nz)
        logger.debug(f"    Grid divisors: X={x_divisors}, Y={y_divisors}, Z={z_divisors}")
        
        best_decomp = None
        best_score = float('inf')
        
        # Try all combinations of divisors that multiply to required_ranks
        for px in rank_divisors:
            if required_ranks % px != 0:
                continue
            
            remaining_y_z = required_ranks // px
            
            for py in rank_divisors:
                if remaining_y_z % py != 0:
                    continue
                
                pz = remaining_y_z // py
                
                # Check if this decomposition is compatible with grid
                x_compatible = (nx % px == 0) if px > 0 else True
                y_compatible = (ny % py == 0) if py > 0 else True
                z_compatible = (nz % pz == 0) if pz > 0 else True
                
                if x_compatible and y_compatible and z_compatible:
                    # Perfect match - return immediately
                    logger.info(f"    Found perfect decomposition: {px} × {py} × {pz}")
                    return (px, py, pz)
                
                # Score based on how "incompatible" it is
                # We'll prefer decompositions that are mostly compatible
                incompatibility = 0
                if not x_compatible:
                    incompatibility += abs(nx % px) * 100
                if not y_compatible:
                    incompatibility += abs(ny % py) * 100
                if not z_compatible:
                    incompatibility += abs(nz % pz) * 100
                
                if incompatibility < best_score:
                    best_score = incompatibility
                    best_decomp = (px, py, pz)
        
        if best_decomp:
            px, py, pz = best_decomp
            if best_score == 0:
                logger.info(f"    Found compatible decomposition: {px} × {py} × {pz}")
            else:
                logger.warning(f"    Using best-effort decomposition {px} × {py} × {pz} "
                              f"(incompatibility score: {best_score})")
            return best_decomp
        
        # Fallback: use the approach from the original code
        if self.initial_procs:
            # Use the initial processor decomposition as a base
            base_x, base_y, base_z = self.initial_procs
            logger.info(f"    Base MPI decomposition: {base_x} × {base_y} × {base_z}")
            
            # Calculate scaling factors to reach required ranks
            # We want to maintain the same ratio as the initial decomposition
            total_base = base_x * base_y * base_z
            if total_base > 0:
                # Scale each dimension proportionally
                scale_factor = (required_ranks / total_base) ** (1/3)
                
                # Try to find integer factors that multiply to required_ranks
                # and are close to the scaled values
                new_x = max(1, round(base_x * scale_factor))
                new_y = max(1, round(base_y * scale_factor))
                new_z = max(1, round(base_z * scale_factor))
                
                # Adjust to ensure we get exactly the required number of ranks
                # We'll adjust one dimension at a time to get closer to required_ranks
                actual_ranks = new_x * new_y * new_z
                
                # Iteratively adjust to get closer to required_ranks
                adjustments = 0
                max_adjustments = 20
                
                while actual_ranks != required_ranks and adjustments < max_adjustments:
                    if actual_ranks < required_ranks:
                        # Need more ranks, increase the largest dimension
                        if new_x >= new_y and new_x >= new_z:
                            new_x += 1
                        elif new_y >= new_x and new_y >= new_z:
                            new_y += 1
                        else:
                            new_z += 1
                    else:
                        # Need fewer ranks, decrease the largest dimension (but not below 1)
                        if new_x >= new_y and new_x >= new_z and new_x > 1:
                            new_x -= 1
                        elif new_y >= new_x and new_y >= new_z and new_y > 1:
                            new_y -= 1
                        elif new_z > 1:
                            new_z -= 1
                        else:
                            # Can't decrease further, break
                            break
                    
                    actual_ranks = new_x * new_y * new_z
                    adjustments += 1
                
                # If we still don't have the exact number, try a different approach
                if actual_ranks != required_ranks:
                    # Try to find factors of required_ranks
                    factors = self._find_factors(required_ranks)
                    if len(factors) >= 3:
                        new_x, new_y, new_z = factors[0], factors[1], factors[2]
                        actual_ranks = new_x * new_y * new_z
                
                logger.info(f"    Adjusted MPI decomposition: {new_x} × {new_y} × {new_z} = {actual_ranks}")
                
                # Update configuration
                cfg['nproc_x'] = new_x
                cfg['nproc_y'] = new_y
                cfg['nproc_z'] = new_z
            else:
                # Fallback: distribute evenly across dimensions
                cube_root = round(required_ranks ** (1/3))
                new_x = new_y = new_z = cube_root
                
                # Adjust to get exactly required_ranks
                actual_ranks = new_x * new_y * new_z
                if actual_ranks != required_ranks:
                    # Simple adjustment - increase one dimension
                    new_x = required_ranks // (new_y * new_z)
                    if new_x * new_y * new_z != required_ranks:
                        # Try another approach
                        new_y = required_ranks // (new_x * new_z)
                        if new_x * new_y * new_z != required_ranks:
                            new_z = required_ranks // (new_x * new_y)
                
                logger.info(f"    Fallback MPI decomposition: {new_x} × {new_y} × {new_z} = {new_x * new_y * new_z}")
                cfg['nproc_x'] = new_x
                cfg['nproc_y'] = new_y
                cfg['nproc_z'] = new_z
        else:
            # No initial processor decomposition provided, distribute evenly
            cube_root = round(required_ranks ** (1/3))
            new_x = new_y = new_z = cube_root
            
            # Adjust to get exactly required_ranks
            actual_ranks = new_x * new_y * new_z
            if actual_ranks != required_ranks:
                # Simple adjustment - increase one dimension
                new_x = required_ranks // (new_y * new_z)
                if new_x * new_y * new_z != required_ranks:
                    # Try another approach
                    new_y = required_ranks // (new_x * new_z)
                    if new_x * new_y * new_z != required_ranks:
                        new_z = required_ranks // (new_x * new_y)
            
            logger.info(f"    Default MPI decomposition: {new_x} × {new_y} × {new_z} = {new_x * new_y * new_z}")
            cfg['nproc_x'] = new_x
            cfg['nproc_y'] = new_y
            cfg['nproc_z'] = new_z
        
        return cfg
    
    def _get_divisors(self, n):
        """
        Get all divisors of n.
        
        Args:
            n: Integer to find divisors for
            
        Returns:
            Sorted list of divisors
        """
        if n <= 0:
            return [1]
        
        divisors = []
        for i in range(1, int(n**0.5) + 1):
            if n % i == 0:
                divisors.append(i)
                if i != n // i:
                    divisors.append(n // i)
        
        return sorted(divisors)
    
    def _find_factors(self, n):
        """
        Find three factors of n that are as close to each other as possible.
        """
        # Start with cube root
        a = round(n ** (1/3))
        
        # Find the closest factor to a
        while a > 0 and n % a != 0:
            a -= 1
        
        if a == 0:
            a = 1
        
        # Now find factors of n//a
        remaining = n // a
        b = round(remaining ** (1/2))
        
        # Find the closest factor to b
        while b > 0 and remaining % b != 0:
            b -= 1
        
        if b == 0:
            b = 1
        
        c = remaining // b
        
        return [a, b, c]
    
    def _scale_to_multiple(self, value, divisor):
        """
        Scale a value to the nearest multiple of divisor.
        Returns the nearest multiple that is >= value.
        
        Examples:
            _scale_to_multiple(840, 224) → 896  (224 * 4 = 896)
            _scale_to_multiple(480, 64)  → 512  (64 * 8 = 512)
        """
        if divisor <= 0:
            return value
        
        # Round up to nearest multiple
        multiples_needed = (value + divisor - 1) // divisor
        return multiples_needed * divisor
    
    def _ensure_divisibility_with_grid_scaling(self, cfg, required_ranks):
        """
        Verify that MPI decomposition divides evenly into grid dimensions.
        If not divisible, scale grid dimensions to make it so.
        """
        nx = safe_int(cfg.get('nx')) or 1
        ny = safe_int(cfg.get('ny')) or 1
        nz = safe_int(cfg.get('nz')) or 1
        
        nproc_x = safe_int(cfg.get('nproc_x')) or 1
        nproc_y = safe_int(cfg.get('nproc_y')) or 1
        nproc_z = safe_int(cfg.get('nproc_z')) or 1
        
        Lx = safe_float(cfg.get('Lx')) or 1.0
        Ly = safe_float(cfg.get('Ly')) or 1.0
        Lz = safe_float(cfg.get('Lz')) or 1.0
        
        # Check and fix X dimension
        if nx % nproc_x != 0:
            nx_old = nx
            nx_new = self._scale_to_multiple(nx, nproc_x)
            Lx_new = Lx * (nx_new / nx_old) if nx_old > 0 else Lx
            cfg['nx'] = nx_new
            cfg['Lx'] = Lx_new
            logger.warning(f"  Grid X scaled for divisibility: {nx_old} → {nx_new} cells, {Lx} → {Lx_new} domain")
        
        # Check and fix Y dimension
        if ny % nproc_y != 0:
            ny_old = ny
            ny_new = self._scale_to_multiple(ny, nproc_y)
            Ly_new = Ly * (ny_new / ny_old) if ny_old > 0 else Ly
            cfg['ny'] = ny_new
            cfg['Ly'] = Ly_new
            logger.warning(f"  Grid Y scaled for divisibility: {ny_old} → {ny_new} cells, {Ly} → {Ly_new} domain")
        
        # Check and fix Z dimension
        if nz % nproc_z != 0:
            nz_old = nz
            nz_new = self._scale_to_multiple(nz, nproc_z)
            Lz_new = Lz * (nz_new / nz_old) if nz_old > 0 else Lz
            cfg['nz'] = nz_new
            cfg['Lz'] = Lz_new
            logger.warning(f"  Grid Z scaled for divisibility: {nz_old} → {nz_new} cells, {Lz} → {Lz_new} domain")
    
    def print_scaling_table(self):
        """Print color-coded verification table with active dimension highlighting."""
        print("\n" + "=" * 140)
        print(f"{Colors.BOLD}{Colors.CYAN}WEAK SCALING VERIFICATION TABLE{Colors.RESET}")
        print(f"Scaling Mode: {Colors.BOLD}{self.dims}D{Colors.RESET}")
        print("=" * 140)
        
        # Header
        header = (
            f"{'Node':>5} │ "
            f"{'MPI':>5} │ "
            f"{'Active':>8} │ "
            f"{Colors.BOLD}Domain (Lx, Ly, Lz){Colors.RESET:>30} │ "
            f"{Colors.BOLD}Cells (nx, ny, nz){Colors.RESET:>30} │ "
            f"{Colors.BOLD}MPI (px, py, pz){Colors.RESET:>30} │ "
            f"{Colors.BOLD}PPC (x, y, z){Colors.RESET:>20}"
        )
        print(header)
        print("─" * 140)
        
        # Data rows
        for node_count, cfg in self.all_configs:
            active = cfg.get('active_dim', 'BASELINE')
            
            # Colorize active dimension
            def colorize(value, dim):
                if active == dim and active != 'BASELINE':
                    return f"{Colors.RED}{value}{Colors.RESET}"
                return str(value)
            
            # Format values
            lx_val = f"{cfg.get('Lx'):.1f}" if isinstance(cfg.get('Lx'), (int, float)) else str(cfg.get('Lx', 'N/A'))
            ly_val = f"{cfg.get('Ly'):.1f}" if isinstance(cfg.get('Ly'), (int, float)) else str(cfg.get('Ly', 'N/A'))
            lz_val = f"{cfg.get('Lz'):.1f}" if isinstance(cfg.get('Lz'), (int, float)) else str(cfg.get('Lz', 'N/A'))
            domain_str = f"({colorize(lx_val, 'X')}, {colorize(ly_val, 'Y')}, {colorize(lz_val, 'Z')})"
            
            nx_val = cfg.get('nx', 'N/A')
            ny_val = cfg.get('ny', 'N/A')
            nz_val = cfg.get('nz', 'N/A')
            cells_str = f"({colorize(nx_val, 'X')}, {colorize(ny_val, 'Y')}, {colorize(nz_val, 'Z')})"
            
            nproc_x_val = cfg.get('nproc_x', 'N/A')
            nproc_y_val = cfg.get('nproc_y', 'N/A')
            nproc_z_val = cfg.get('nproc_z', 'N/A')
            mpi_str = f"({colorize(nproc_x_val, 'X')}, {colorize(nproc_y_val, 'Y')}, {colorize(nproc_z_val, 'Z')})"
            
            # Particles per cell (should be constant)
            ppc_x = cfg.get('num_particles_x', 'N/A')
            ppc_y = cfg.get('num_particles_y', 'N/A')
            ppc_z = cfg.get('num_particles_z', 'N/A')
            ppc_str = f"({ppc_x}, {ppc_y}, {ppc_z})"
            
            # Active dimension display
            if active == 'BASELINE':
                active_display = f"{Colors.GREEN}{active}{Colors.RESET}"
            else:
                active_display = f"{Colors.RED}{Colors.BOLD}{active}{Colors.RESET}"
            
            # Handle None values and non-numeric values for MPI process counts
            nproc_x = safe_int(cfg.get('nproc_x')) or 1
            nproc_y = safe_int(cfg.get('nproc_y')) or 1
            nproc_z = safe_int(cfg.get('nproc_z')) or 1
            total_ranks = nproc_x * nproc_y * nproc_z
            
            print(
                f"{node_count:>5} │ "
                f"{total_ranks:>5} │ "
                f"{active_display:>8} │ "
                f"{domain_str:>30} │ "
                f"{cells_str:>30} │ "
                f"{mpi_str:>30} │ "
                f"{ppc_str:>20}"
            )
        
        print("=" * 140)
        
        # Summary
        self._print_summary()
    
    def _print_summary(self):
        """Print summary statistics."""
        print(f"\n{Colors.BOLD}Scaling Pattern Summary:{Colors.RESET}")
        
        if self.dims == 1:
            print(f"  • 1D Scaling: {Colors.RED}X only{Colors.RESET}, Y and Z constant")
        elif self.dims == 2:
            print(f"  • 2D Scaling: {Colors.RED}X→Y→X→Y{Colors.RESET} alternating, Z constant")
        else:
            print(f"  • 3D Scaling: {Colors.RED}X→Y→Z→X→Y→Z{Colors.RESET} cycling")
        
        # Check if particles per cell exist and are constant
        if self.all_configs and 'num_particles_x' in self.all_configs[0][1]:
            base_cfg = self.all_configs[0][1]
            all_constant = all(
                cfg.get('num_particles_x') == base_cfg.get('num_particles_x') and
                cfg.get('num_particles_y') == base_cfg.get('num_particles_y') and
                cfg.get('num_particles_z') == base_cfg.get('num_particles_z')
                for _, cfg in self.all_configs
            )
            
            if all_constant:
                print(f"  • Particles per cell: {Colors.GREEN}✓ CONSTANT{Colors.RESET} "
                      f"({base_cfg.get('num_particles_x')}, {base_cfg.get('num_particles_y')}, {base_cfg.get('num_particles_z')})")
            else:
                print(f"  • Particles per cell: {Colors.RED}✗ NOT CONSTANT{Colors.RESET}")
        
        print()
    
    def generate_job_configs(self):
        """Generate job configurations for all nodes."""
        from core.config import JobConfig
        
        # Generate all cases first
        self.generate_all_cases()
        
        # Convert configs to JobConfig objects
        job_configs = []
        for node_count, cfg in self.all_configs:
            # Create job ID
            job_id = f"node{node_count}"
            
            # Calculate total procs
            num_procs = cfg.get('nproc_x', 1) * cfg.get('nproc_y', 1) * cfg.get('nproc_z', 1)
            
            # Create procs decomposition
            procs_decomposition = (cfg.get('nproc_x', 1), cfg.get('nproc_y', 1), cfg.get('nproc_z', 1))
            
            # Create domain size
            domain_size = (cfg.get('Lx', 1.0), cfg.get('Ly', 1.0), cfg.get('Lz', 1.0))
            
            # Create cell count
            cell_count = (cfg.get('nx', 1), cfg.get('ny', 1), cfg.get('nz', 1))
            
            # Create JobConfig
            job_config = JobConfig(
                job_id=job_id,
                num_nodes=node_count,
                num_procs=num_procs,
                procs_decomposition=procs_decomposition,
                domain_size=domain_size,
                cell_count=cell_count
            )
            
            job_configs.append(job_config)
        
        return job_configs

# ------------------- CLI -------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="1D/2D/3D Weak Scaling Engine with Color-Coded Table",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i config.yaml -o output/ -n 1 2 4 8 -p 112 -d 2
  %(prog)s -i baseline.yml -o scaled/ -n 1 2 4 8 16 -p 128 -k 2 -d 3
        """
    )
    parser.add_argument('-i', '--input', required=True, 
                        help="Input YAML file (.yaml or .yml)")
    parser.add_argument('-o', '--output', required=True, 
                        help="Output directory for generated files")
    parser.add_argument('-n', '--nodes', nargs='+', type=int, required=True, 
                        help="Node counts (e.g., 1 2 4 8)")
    parser.add_argument('-p', '--procs_per_node', type=int, default=1, 
                        help="Processes per node (default: 1)")
    parser.add_argument('-k', '--scale_factor', type=float, default=2.0, 
                        help="Weak scaling factor (default: 2.0)")
    parser.add_argument('-d', '--dims', type=int, default=3, choices=[1, 2, 3], 
                        help="Scaling dimensions: 1=X only, 2=X→Y, 3=X→Y→Z (default: 3)")
    args = parser.parse_args()

    # Validate YAML extension
    if not args.input.lower().endswith(('.yaml', '.yml')):
        logger.error(f"Input file must be YAML (.yaml or .yml), got: {args.input}")
        raise ValueError(f"Invalid file extension: {args.input}")
    
    # Check if input file exists
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        raise FileNotFoundError(args.input)
    
    logger.info(f"{Colors.BOLD}{Colors.CYAN}Starting Weak Scaling Engine{Colors.RESET}")
    logger.info("=" * 80)

    # Load YAML config
    with open(args.input, 'r') as f:
        yaml_config = yaml.safe_load(f)
    
    # Extract parameters from YAML
    initial_domain = yaml_config.get('initial_domain')
    initial_cells = yaml_config.get('initial_cells')
    initial_procs = yaml_config.get('initial_procs')
    particles_per_cell = None
    if 'particles_per_cell' in yaml_config:
        ppc = yaml_config['particles_per_cell']
        particles_per_cell = [ppc.get('x'), ppc.get('y'), ppc.get('z')]
    
    engine = ScalingEngine(
        input_file=yaml_config.get('input_file', 'inputfiles/os-stdin'),
        output_dir=args.output,
        nodes=args.nodes,
        procs_per_node=args.procs_per_node,
        scale_factor=args.scale_factor,
        dims=args.dims,
        initial_domain=initial_domain,
        initial_cells=initial_cells,
        initial_procs=initial_procs,
        particles_per_cell=particles_per_cell
    )
    
    try:
        engine.generate_all_cases()
        logger.info(f"{Colors.GREEN}{Colors.BOLD}Scaling completed successfully!{Colors.RESET}")
    except Exception as e:
        logger.error(f"{Colors.RED}Scaling failed: {e}{Colors.RESET}")
        raise
