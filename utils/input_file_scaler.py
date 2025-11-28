#!/usr/bin/env python3
"""
Standalone input file scaler for PIC simulations.
This tool generates scaled input files for weak scaling tests.
"""

import os
import re
import yaml
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

class InputFileScaler:
    """Scale PIC input files for weak scaling tests."""
    
    def __init__(self, template_file: str, output_dir: str):
        self.template_file = Path(template_file)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Read template content
        with open(self.template_file, 'r') as f:
            self.template_content = f.read()
    
    def generate_scaled_file(
        self,
        node_count: int,
        domain_size: Tuple[float, float, float],
        cell_count: Tuple[int, int, int],
        mpi_decomp: Tuple[int, int, int],
        particles_per_cell: Tuple[int, int, int]
    ) -> str:
        """
        Generate a scaled input file for a specific node configuration.
        
        Args:
            node_count: Node number
            domain_size: (Lx, Ly, Lz) values
            cell_count: (nx, ny, nz) values
            mpi_decomp: (nproc_x, nproc_y, nproc_z) values
            particles_per_cell: (npcelx, npcely, npcelz) values
            
        Returns:
            Path to generated file
        """
        # Start with template content
        content = self.template_content
        
        # Extract variable names from content
        Lx = domain_size[0]
        Ly = domain_size[1]
        Lz = domain_size[2]
        
        nx = cell_count[0]
        ny = cell_count[1]
        nz = cell_count[2]
        
        nproc_x = mpi_decomp[0]
        nproc_y = mpi_decomp[1]
        nproc_z = mpi_decomp[2]
        
        npcelx = particles_per_cell[0]
        npcely = particles_per_cell[1]
        npcelz = particles_per_cell[2]
        
        # Replace domain size variables
        content = self._replace_variable(content, 'Lx', Lx)
        content = self._replace_variable(content, 'Ly', Ly)
        content = self._replace_variable(content, 'Lz', Lz)
        
        # Replace cell count variables (note: nx->nxc, etc.)
        content = self._replace_variable(content, 'nxc', nx)
        content = self._replace_variable(content, 'nyc', ny)
        content = self._replace_variable(content, 'nzc', nz)
        
        # Replace MPI decomposition variables
        content = self._replace_variable(content, 'XLEN', nproc_x)
        content = self._replace_variable(content, 'YLEN', nproc_y)
        content = self._replace_variable(content, 'ZLEN', nproc_z)
        
        # Replace particles per cell variables
        content = self._replace_variable(content, 'npcelx', npcelx)
        content = self._replace_variable(content, 'npcely', npcely)
        content = self._replace_variable(content, 'npcelz', npcelz)
        
        # Write output file
        output_filename = f"{self.template_file.stem}_node{node_count}.in"
        output_path = self.output_dir / output_filename
        
        with open(output_path, 'w') as f:
            f.write(content)
        
        logger.info(f"Generated scaled input file: {output_path}")
        return str(output_path)
    
    def _replace_variable(self, content: str, var_name: str, new_value: Any) -> str:
        """
        Replace variable value in content.
        
        Args:
            content: File content
            var_name: Variable name
            new_value: New value
            
        Returns:
            Modified content
        """
        # Build replacement patterns - match ANY value (not just numeric)
        # Use word boundary or non-word character to avoid partial matches
        patterns = [
            # Pattern 1: var = value or var: value
            (rf'({re.escape(var_name)}\s*[=:])\s*\S+', rf'\g<1> {new_value}'),
            # Pattern 2: "var" = value or 'var': value  
            (rf'(["\']?{re.escape(var_name)}["\']?\s*[=:])\s*\S+', rf'\g<1> {new_value}'),
            # Pattern 3: <var> value
            (rf'(<{re.escape(var_name)}>)\s+\S+', rf'\g<1> {new_value}'),
        ]
        
        for pattern, replacement in patterns:
            new_content = re.sub(pattern, replacement, content, flags=re.IGNORECASE | re.MULTILINE, count=1)
            if new_content != content:
                return new_content
        
        # If no pattern matched, log warning but don't fail
        logger.debug(f"Could not find pattern to replace '{var_name}' in input file")
        return content

def main():
    parser = argparse.ArgumentParser(
        description="Scale PIC input files for weak scaling tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --template inputfiles/os-stdin --output output/ --config run.yaml
        """
    )
    
    parser.add_argument(
        "--template",
        required=True,
        help="Template input file to scale"
    )
    
    parser.add_argument(
        "--output",
        default="./scaled_inputs",
        help="Output directory for generated files"
    )
    
    parser.add_argument(
        "--config",
        required=True,
        help="YAML configuration file with scaling parameters"
    )
    
    args = parser.parse_args()
    
    # Load configuration
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)
    
    # Extract scaling parameters
    nodes = config.get('nodes', [1, 2, 4, 8])
    dims = config.get('scaling_dimensions', 2)
    scale_factor = config.get('scaling_factor', 2.0)
    procs_per_node = config.get('procs_per_node', 112)
    
    # Extract initial values
    initial_domain = config.get('initial_domain', [84.0, 48.0, 1.0])
    initial_cells = config.get('initial_cells', [840, 480, 1])
    initial_procs = config.get('initial_procs', [14, 8, 1])
    particles_per_cell = config.get('particles_per_cell', [20, 20, 1])
    
    # Create scaler
    scaler = InputFileScaler(args.template, args.output)
    
    # Generate files for each node
    print(f"\n{Colors.BOLD}{Colors.CYAN}GENERATING SCALED INPUT FILES{Colors.RESET}")
    print("=" * 80)
    
    # Store configurations for table
    configs = []
    
    for i, node_count in enumerate(nodes):
        if i == 0:
            # Node 1 - baseline
            Lx, Ly, Lz = initial_domain
            nx, ny, nz = initial_cells
            nproc_x, nproc_y, nproc_z = initial_procs
            npcelx, npcely, npcelz = particles_per_cell
            active_dim = 'BASELINE'
        else:
            # Determine active dimension
            if dims == 1:
                active_dim = 'X'
            elif dims == 2:
                active_dim = 'X' if i % 2 == 1 else 'Y'
            else:  # dims == 3
                cycle = ['X', 'Y', 'Z']
                active_dim = cycle[(i - 1) % 3]
            
            # Start with previous values
            prev_config = configs[-1]
            Lx, Ly, Lz = prev_config['domain']
            nx, ny, nz = prev_config['cells']
            nproc_x, nproc_y, nproc_z = prev_config['mpi']
            npcelx, npcely, npcelz = particles_per_cell
            
            # Scale active dimension
            if active_dim == 'X':
                Lx *= scale_factor
                nx = int(nx * scale_factor)
                nproc_x = int(nproc_x * scale_factor)
            elif active_dim == 'Y':
                Ly *= scale_factor
                ny = int(ny * scale_factor)
                nproc_y = int(nproc_y * scale_factor)
            elif active_dim == 'Z':
                Lz *= scale_factor
                nz = int(nz * scale_factor)
                nproc_z = int(nproc_z * scale_factor)
            
            # MPI correction
            required_ranks = node_count * procs_per_node
            actual_ranks = nproc_x * nproc_y * nproc_z
            
            if actual_ranks != required_ranks:
                if active_dim == 'X':
                    nproc_x = required_ranks // (nproc_y * nproc_z)
                elif active_dim == 'Y':
                    nproc_y = required_ranks // (nproc_x * nproc_z)
                elif active_dim == 'Z':
                    nproc_z = required_ranks // (nproc_x * nproc_y)
        
        # Store configuration
        config_entry = {
            'node': node_count,
            'active': active_dim,
            'domain': (Lx, Ly, Lz),
            'cells': (nx, ny, nz),
            'mpi': (nproc_x, nproc_y, nproc_z),
            'particles': (npcelx, npcely, npcelz)
        }
        configs.append(config_entry)
        
        # Generate file
        scaler.generate_scaled_file(
            node_count,
            (Lx, Ly, Lz),
            (nx, ny, nz),
            (nproc_x, nproc_y, nproc_z),
            (npcelx, npcely, npcelz)
        )
    
    # Print verification table
    print("\n" + "=" * 120)
    print(f"{Colors.BOLD}{Colors.CYAN}WEAK SCALING VERIFICATION TABLE{Colors.RESET}")
    print(f"Scaling Mode: {Colors.BOLD}{dims}D{Colors.RESET}")
    print("=" * 120)
    
    # Header
    header = (
        f"{'Node':>5} │ "
        f"{'MPI':>5} │ "
        f"{'Active':>8} │ "
        f"{Colors.BOLD}Domain (Lx, Ly, Lz){Colors.RESET:>25} │ "
        f"{Colors.BOLD}Cells (nx, ny, nz){Colors.RESET:>25} │ "
        f"{Colors.BOLD}MPI (px, py, pz){Colors.RESET:>20}"
    )
    print(header)
    print("─" * 120)
    
    # Data rows
    for config in configs:
        node_count = config['node']
        active = config['active']
        Lx, Ly, Lz = config['domain']
        nx, ny, nz = config['cells']
        nproc_x, nproc_y, nproc_z = config['mpi']
        
        # Colorize active dimension
        def colorize(value, dim):
            if active == dim and active != 'BASELINE':
                return f"{Colors.RED}{value}{Colors.RESET}"
            return str(value)
        
        # Format values
        domain_str = f"({colorize(f'{Lx:.0f}', 'X')}, {colorize(f'{Ly:.0f}', 'Y')}, {colorize(f'{Lz:.0f}', 'Z')})"
        cells_str = f"({colorize(nx, 'X')}, {colorize(ny, 'Y')}, {colorize(nz, 'Z')})"
        mpi_str = f"({colorize(nproc_x, 'X')}, {colorize(nproc_y, 'Y')}, {colorize(nproc_z, 'Z')})"
        
        # Active dimension display
        if active == 'BASELINE':
            active_display = f"{Colors.GREEN}{active}{Colors.RESET}"
        else:
            active_display = f"{Colors.RED}{Colors.BOLD}{active}{Colors.RESET}"
        
        total_ranks = nproc_x * nproc_y * nproc_z
        
        print(
            f"{node_count:>5} │ "
            f"{total_ranks:>5} │ "
            f"{active_display:>8} │ "
            f"{domain_str:>25} │ "
            f"{cells_str:>25} │ "
            f"{mpi_str:>20}"
        )
    
    print("─" * 120)
    print(f"\n{Colors.BOLD}Scaling Pattern Summary:{Colors.RESET}")
    if dims == 1:
        pattern = "1D Scaling: X only"
    elif dims == 2:
        pattern = "2D Scaling: X→Y→X→Y alternating, Z constant"
    else:
        pattern = "3D Scaling: X→Y→Z cycling"
    print(f"  • {pattern}")
    print(f"  • Particles per cell: {Colors.GREEN}✓ CONSTANT{Colors.RESET} ({particles_per_cell[0]}, {particles_per_cell[1]}, {particles_per_cell[2]})")
    print("=" * 120)

if __name__ == "__main__":
    main()