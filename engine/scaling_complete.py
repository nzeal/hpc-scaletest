#!/usr/bin/env python3
"""
scaling.py - General 1D/2D/3D weak scaling engine for PIC simulation inputs
with color-coded table highlighting and detailed logging.

Features:
- Node 1 uses baseline parameters from YAML
- Weak scaling applied step-by-step per dimension
- MPI decomposition corrected to match total ranks
- Particles-per-cell remain constant
- Supports arbitrary YAML input via -i flag
- YAML extension validation (.yaml or .yml)
- Color-coded table output with active dimension highlighting
- Detailed step-by-step logging
"""

import os
import yaml
import copy
import argparse
import logging
from pathlib import Path
from typing import Dict, Any, Tuple, List, Optional

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
    return value

def scale_int(value, factor):
    num = safe_int(value)
    if num is not None:
        return int(round(num * factor))
    return value

# ------------------- Scaling Engine -------------------
class ScalingEngine:
    def __init__(self, input_file, output_dir, nodes, procs_per_node, scale_factor=2.0, dims=3, 
                 initial_domain=None, initial_cells=None, initial_procs=None, particles_per_cell=None):
        logger.info(f"Initializing Scaling Engine")
        logger.info(f"  Input file: {input_file}")
        logger.info(f"  Output directory: {output_dir}")
        logger.info(f"  Node sequence: {nodes}")
        logger.info(f"  Procs per node: {procs_per_node}")
        logger.info(f"  Scale factor: {scale_factor}")
        logger.info(f"  Scaling dimensions: {dims}D")
        
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

        # Load baseline config - handle both YAML and plain text input files
        logger.info(f"Loading baseline configuration from {input_file}")
        if input_file.lower().endswith(('.yaml', '.yml')):
            # YAML file - load as before
            with open(input_file, 'r') as f:
                self.base_config = yaml.safe_load(f)
            logger.info(f"  Loaded YAML parameters: {list(self.base_config.keys())}")
        else:
            # Plain text input file - use input analyzer
            from utils.input_analyzer import InputFileAnalyzer
            analyzer = InputFileAnalyzer(Path(input_file))
            params = analyzer.params
            
            # Convert to the format expected by scaling engine
            self.base_config = {
                'Lx': params.Lx,
                'Ly': params.Ly,
                'Lz': params.Lz,
                'nx': params.nxc,
                'ny': params.nyc,
                'nz': params.nzc,
                'nproc_x': params.XLEN,
                'nproc_y': params.YLEN,
                'nproc_z': params.ZLEN,
                'num_particles_x': params.npcelx,
                'num_particles_y': params.npcely,
                'num_particles_z': params.npcelz
            }
            
            # Remove None values
            self.base_config = {k: v for k, v in self.base_config.items() if v is not None}
            logger.info(f"  Loaded text parameters: {list(self.base_config.keys())}")
        
        # Override with initial values if provided
        if initial_domain:
            self.base_config['Lx'] = initial_domain[0]
            self.base_config['Ly'] = initial_domain[1]
            self.base_config['Lz'] = initial_domain[2]
            logger.info(f"  Overriding domain with initial values: {initial_domain}")
        
        if initial_cells:
            self.base_config['nx'] = initial_cells[0]
            self.base_config['ny'] = initial_cells[1]
            self.base_config['nz'] = initial_cells[2]
            logger.info(f"  Overriding cells with initial values: {initial_cells}")
        
        if initial_procs:
            self.base_config['nproc_x'] = initial_procs[0]
            self.base_config['nproc_y'] = initial_procs[1]
            self.base_config['nproc_z'] = initial_procs[2]
            logger.info(f"  Overriding MPI decomposition with initial values: {initial_procs}")
            
        # Set particles per cell if provided
        if particles_per_cell:
            self.base_config['num_particles_x'] = particles_per_cell[0]
            self.base_config['num_particles_y'] = particles_per_cell[1]
            self.base_config['num_particles_z'] = particles_per_cell[2]
            logger.info(f"  Setting particles per cell: {particles_per_cell}")

        self.output_basename = os.path.splitext(os.path.basename(input_file))[0]
        os.makedirs(output_dir, exist_ok=True)
        
        # Store all node configs for table
        self.all_configs = []

    def generate_all_cases(self):
        logger.info("\nGenerating scaling configurations for all nodes...")
        logger.info("=" * 80)
        
        for idx, node_count in enumerate(self.nodes):
            logger.info(f"\nProcessing Step {idx}: Node {node_count}")
            
            cfg = self.scale_node(idx, node_count)
            self.all_configs.append((node_count, cfg))
            
            out_file = os.path.join(self.output_dir, f"{self.output_basename}_node{node_count}.yaml")
            with open(out_file, 'w') as f:
                yaml.dump(cfg, f)
            logger.info(f"  ✓ Generated: {out_file}")
        
        logger.info("\n" + "=" * 80)
        logger.info("All configurations generated successfully\n")
        
        # Print color-coded verification table
        self.print_scaling_table()

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

    def scale_node(self, step_index, node_count):
        cfg = copy.deepcopy(self.base_config)
        
        if step_index == 0:
            # Node 1 baseline, no scaling
            logger.info(f"  {Colors.GREEN}BASELINE{Colors.RESET} - Using exact values from YAML")
            logger.info(f"    Domain: Lx={cfg.get('Lx', 'N/A')}, Ly={cfg.get('Ly', 'N/A')}, Lz={cfg.get('Lz', 'N/A')}")
            logger.info(f"    Cells:  nx={cfg.get('nx', 'N/A')}, ny={cfg.get('ny', 'N/A')}, nz={cfg.get('nz', 'N/A')}")
            logger.info(f"    MPI:    px={cfg.get('nproc_x', 'N/A')}, py={cfg.get('nproc_y', 'N/A')}, pz={cfg.get('nproc_z', 'N/A')}")
            cfg['active_dim'] = 'BASELINE'
            return cfg

        # Determine active dimension
        active_dim = self.get_active_dimension(step_index)
        factor = self.scale_factor
        
        logger.info(f"  {Colors.RED}Scaling dimension: {active_dim}{Colors.RESET} (factor={factor})")

        # Log before scaling
        before = {
            'Lx': cfg.get('Lx'), 'Ly': cfg.get('Ly'), 'Lz': cfg.get('Lz'),
            'nx': cfg.get('nx'), 'ny': cfg.get('ny'), 'nz': cfg.get('nz'),
            'nproc_x': cfg.get('nproc_x'), 'nproc_y': cfg.get('nproc_y'), 'nproc_z': cfg.get('nproc_z')
        }

        # Apply scaling to active dimension only
        # Only scale if the parameter exists
        if active_dim == 'X':
            if 'Lx' in cfg:
                cfg['Lx'] = scale_value(cfg['Lx'], factor)
            if 'nx' in cfg:
                cfg['nx'] = scale_int(cfg['nx'], factor)
            if 'nproc_x' in cfg:
                cfg['nproc_x'] = scale_int(cfg['nproc_x'], factor)
        elif active_dim == 'Y':
            if 'Ly' in cfg:
                cfg['Ly'] = scale_value(cfg['Ly'], factor)
            if 'ny' in cfg:
                cfg['ny'] = scale_int(cfg['ny'], factor)
            if 'nproc_y' in cfg:
                cfg['nproc_y'] = scale_int(cfg['nproc_y'], factor)
        elif active_dim == 'Z':
            if 'Lz' in cfg:
                cfg['Lz'] = scale_value(cfg['Lz'], factor)
            if 'nz' in cfg:
                cfg['nz'] = scale_int(cfg['nz'], factor)
            if 'nproc_z' in cfg:
                cfg['nproc_z'] = scale_int(cfg['nproc_z'], factor)
        
        # Log after scaling
        logger.info(f"    Domain: Lx={before.get('Lx', 'N/A')}→{cfg.get('Lx', 'N/A')}, Ly={before.get('Ly', 'N/A')}→{cfg.get('Ly', 'N/A')}, Lz={before.get('Lz', 'N/A')}→{cfg.get('Lz', 'N/A')}")
        logger.info(f"    Cells:  nx={before.get('nx', 'N/A')}→{cfg.get('nx', 'N/A')}, ny={before.get('ny', 'N/A')}→{cfg.get('ny', 'N/A')}, nz={before.get('nz', 'N/A')}→{cfg.get('nz', 'N/A')}")
        logger.info(f"    MPI:    px={before.get('nproc_x', 'N/A')}→{cfg.get('nproc_x', 'N/A')}, py={before.get('nproc_y', 'N/A')}→{cfg.get('nproc_y', 'N/A')}, pz={before.get('nproc_z', 'N/A')}→{cfg.get('nproc_z', 'N/A')}")

        # MPI decomposition correction
        cfg = self.correct_mpi(cfg, active_dim, node_count)
        cfg['active_dim'] = active_dim
        
        return cfg

    def get_active_dimension(self, step_index):
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
        required_ranks = node_count * self.procs_per_node
        npx = cfg.get('nproc_x', 1)
        npy = cfg.get('nproc_y', 1)
        npz = cfg.get('nproc_z', 1)

        actual_ranks = npx * npy * npz

        if actual_ranks == required_ranks:
            return cfg  # nothing to do
        
        logger.info(f"  {Colors.YELLOW}MPI correction needed{Colors.RESET}: {actual_ranks} → {required_ranks} ranks")

        # Correct active dimension only
        if active_dim == 'X':
            if 'nproc_x' in cfg:
                npx_new = required_ranks // (npy * npz)
                if npx_new * npy * npz != required_ranks:
                    raise ValueError(f"MPI correction failed for node {node_count}")
                logger.info(f"    Adjusting nproc_x: {npx} → {npx_new}")
                cfg['nproc_x'] = npx_new
        elif active_dim == 'Y':
            if 'nproc_y' in cfg:
                npy_new = required_ranks // (npx * npz)
                if npx * npy_new * npz != required_ranks:
                    raise ValueError(f"MPI correction failed for node {node_count}")
                logger.info(f"    Adjusting nproc_y: {npy} → {npy_new}")
                cfg['nproc_y'] = npy_new
        elif active_dim == 'Z':
            if 'nproc_z' in cfg:
                npz_new = required_ranks // (npx * npy)
                if npx * npy * npz_new != required_ranks:
                    raise ValueError(f"MPI correction failed for node {node_count}")
                logger.info(f"    Adjusting nproc_z: {npz} → {npz_new}")
                cfg['nproc_z'] = npz_new
        return cfg
    
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
            f"{Colors.BOLD}MPI (px, py, pz){Colors.RESET:>30}"
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
            
            # Active dimension display
            if active == 'BASELINE':
                active_display = f"{Colors.GREEN}{active}{Colors.RESET}"
            else:
                active_display = f"{Colors.RED}{Colors.BOLD}{active}{Colors.RESET}"
            
            # Handle None values for MPI process counts
            nproc_x = cfg.get('nproc_x') or 1
            nproc_y = cfg.get('nproc_y') or 1
            nproc_z = cfg.get('nproc_z') or 1
            total_ranks = nproc_x * nproc_y * nproc_z
            
            print(
                f"{node_count:>5} │ "
                f"{total_ranks:>5} │ "
                f"{active_display:>8} │ "
                f"{domain_str:>30} │ "
                f"{cells_str:>30} │ "
                f"{mpi_str:>30}"
            )
        
        print("─" * 140)
        
        # Scaling pattern summary
        if self.dims == 1:
            pattern = "1D Scaling: X only"
        elif self.dims == 2:
            pattern = "2D Scaling: X→Y→X→Y alternating, Z constant"
        else:
            pattern = "3D Scaling: X→Y→Z cycling"
            
        print(f"\n{Colors.BOLD}Scaling Pattern Summary:{Colors.RESET}")
        print(f"  • {pattern}")
        
        # Particles per cell
        if 'num_particles_x' in self.base_config:
            ppc_x = self.base_config['num_particles_x']
            ppc_y = self.base_config.get('num_particles_y', 'N/A')
            ppc_z = self.base_config.get('num_particles_z', 'N/A')
            print(f"  • Particles per cell: {Colors.GREEN}✓ CONSTANT{Colors.RESET} ({ppc_x}, {ppc_y}, {ppc_z})")
        
        print("=" * 140)

# ------------------- CLI Entry Point -------------------
def main():
    parser = argparse.ArgumentParser(
        description="General 1D/2D/3D weak scaling engine for PIC simulation inputs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i run.yaml
  %(prog)s -i config.yaml -o output/
  %(prog)s -i run.yaml --nodes 1 2 4 8
        """
    )
    
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input YAML configuration file (.yaml or .yml)"
    )
    
    parser.add_argument(
        "-o", "--output",
        default="./scaling_configs",
        help="Output directory for generated configurations"
    )
    
    parser.add_argument(
        "--nodes",
        nargs="+",
        type=int,
        default=[1, 2, 4, 8],
        help="Node sequence for scaling (default: 1 2 4 8)"
    )
    
    parser.add_argument(
        "--procs-per-node",
        type=int,
        default=112,
        help="Processors per node (default: 112)"
    )
    
    parser.add_argument(
        "--scale-factor",
        type=float,
        default=2.0,
        help="Scaling factor (default: 2.0)"
    )
    
    parser.add_argument(
        "--dims",
        type=int,
        choices=[1, 2, 3],
        default=2,
        help="Scaling dimensions: 1, 2, or 3 (default: 2)"
    )
    
    args = parser.parse_args()
    
    # Validate input file extension
    if not args.input.lower().endswith(('.yaml', '.yml')):
        raise ValueError("Input file must have .yaml or .yml extension")
    
    # Load configuration
    with open(args.input, 'r') as f:
        config = yaml.safe_load(f)
    
    # Extract initial values from config
    initial_domain = config.get('initial_domain')
    initial_cells = config.get('initial_cells')
    initial_procs = config.get('initial_procs')
    particles_per_cell = config.get('particles_per_cell')
    
    # Create scaling engine
    engine = ScalingEngine(
        input_file=args.input,
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
    
    # Generate all cases
    engine.generate_all_cases()

if __name__ == "__main__":
    main()