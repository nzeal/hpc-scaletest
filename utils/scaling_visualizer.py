#!/usr/bin/env python3
"""
Scaling Visualizer - Color-coded table output for weak scaling verification.

Provides visual feedback on scaling progression with active dimension highlighting.
"""

import logging
from typing import List, Tuple, Any
from dataclasses import dataclass

logger = logging.getLogger(__name__)


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
    
    @staticmethod
    def disable():
        """Disable colors (for non-TTY environments)."""
        Colors.RED = ''
        Colors.GREEN = ''
        Colors.YELLOW = ''
        Colors.BLUE = ''
        Colors.MAGENTA = ''
        Colors.CYAN = ''
        Colors.BOLD = ''
        Colors.UNDERLINE = ''
        Colors.RESET = ''


@dataclass
class ScalingSnapshot:
    """Snapshot of scaling configuration for a single node."""
    node_count: int
    num_procs: int
    active_dimension: str  # 'X', 'Y', 'Z', or 'BASELINE'
    
    # Domain
    domain_x: float
    domain_y: float
    domain_z: float
    
    # Cells
    cells_x: int
    cells_y: int
    cells_z: int
    
    # MPI decomposition
    mpi_x: int
    mpi_y: int
    mpi_z: int
    
    # Particles per cell (constant)
    particles_x: int
    particles_y: int
    particles_z: int


class ScalingVisualizer:
    """Visualize scaling progression with color-coded tables."""
    
    def __init__(self, disable_colors: bool = False):
        """
        Initialize visualizer.
        
        Args:
            disable_colors: If True, disable ANSI colors for non-TTY environments
        """
        if disable_colors:
            Colors.disable()
    
    def print_scaling_table(self, snapshots: List[ScalingSnapshot], scaling_dimensions: int):
        """
        Print color-coded scaling verification table.
        
        Args:
            snapshots: List of scaling snapshots for each node
            scaling_dimensions: Number of scaling dimensions (1, 2, or 3)
        """
        print("\n" + "=" * 140)
        print(f"{Colors.BOLD}{Colors.CYAN}WEAK SCALING VERIFICATION TABLE{Colors.RESET}")
        print(f"Scaling Mode: {Colors.BOLD}{scaling_dimensions}D{Colors.RESET}")
        print("=" * 140)
        
        # Header
        header = (
            f"{'Node':>5} │ "
            f"{'MPI':>5} │ "
            f"{'Active':>8} │ "
            f"{Colors.BOLD}Domain (Lx, Ly, Lz){Colors.RESET:<35} │ "
            f"{Colors.BOLD}Cells (nx, ny, nz){Colors.RESET:<35} │ "
            f"{Colors.BOLD}MPI (px, py, pz){Colors.RESET:<35}"
        )
        print(header)
        print("─" * 140)
        
        # Data rows
        for snapshot in snapshots:
            self._print_row(snapshot)
        
        print("=" * 140)
        
        # Summary
        self._print_summary(snapshots, scaling_dimensions)
    
    def _print_row(self, snap: ScalingSnapshot):
        """Print a single row with color highlighting."""
        active = snap.active_dimension
        
        # Color the active dimension
        def colorize(value, dim):
            if active == dim and active != 'BASELINE':
                return f"{Colors.RED}{value}{Colors.RESET}"
            return str(value)
        
        # Format values
        lx_val = f'{snap.domain_x:.1f}' if isinstance(snap.domain_x, (int, float)) else str(snap.domain_x)
        ly_val = f'{snap.domain_y:.1f}' if isinstance(snap.domain_y, (int, float)) else str(snap.domain_y)
        lz_val = f'{snap.domain_z:.1f}' if isinstance(snap.domain_z, (int, float)) else str(snap.domain_z)
        domain_str = (
            f"({colorize(lx_val, 'X'):>6}, "
            f"{colorize(ly_val, 'Y'):>6}, "
            f"{colorize(lz_val, 'Z'):>6})"
        )
        
        cells_str = (
            f"({colorize(snap.cells_x, 'X'):>6}, "
            f"{colorize(snap.cells_y, 'Y'):>6}, "
            f"{colorize(snap.cells_z, 'Z'):>6})"
        )
        
        mpi_str = (
            f"({colorize(snap.mpi_x, 'X'):>6}, "
            f"{colorize(snap.mpi_y, 'Y'):>6}, "
            f"{colorize(snap.mpi_z, 'Z'):>6})"
        )
        
        # Active dimension display
        if active == 'BASELINE':
            active_display = f"{Colors.GREEN}{active}{Colors.RESET}"
        else:
            active_display = f"{Colors.RED}{Colors.BOLD}{active}{Colors.RESET}"
        
        print(
            f"{snap.node_count:>5} │ "
            f"{snap.num_procs:>5} │ "
            f"{active_display:>8} │ "
            f"{domain_str:<35} │ "
            f"{cells_str:<35} │ "
            f"{mpi_str:<35}"
        )
    
    def _print_summary(self, snapshots: List[ScalingSnapshot], dims: int):
        """Print summary statistics."""
        print()
        print(f"{Colors.BOLD}Scaling Pattern Summary:{Colors.RESET}")
        
        if dims == 1:
            print(f"  • 1D Scaling: {Colors.RED}X only{Colors.RESET}, Y and Z constant")
        elif dims == 2:
            print(f"  • 2D Scaling: {Colors.RED}X→Y→X→Y{Colors.RESET} alternating, Z constant")
        else:
            print(f"  • 3D Scaling: {Colors.RED}X→Y→Z→X→Y→Z{Colors.RESET} cycling")
        
        # Verify particles per cell constant
        base = snapshots[0]
        all_constant = all(
            s.particles_x == base.particles_x and
            s.particles_y == base.particles_y and
            s.particles_z == base.particles_z
            for s in snapshots
        )
        
        if all_constant:
            print(f"  • Particles per cell: {Colors.GREEN}✓ CONSTANT{Colors.RESET} "
                  f"({base.particles_x}, {base.particles_y}, {base.particles_z})")
        else:
            print(f"  • Particles per cell: {Colors.RED}✗ NOT CONSTANT{Colors.RESET}")
        
        # Scaling factors
        if len(snapshots) > 1:
            print(f"\n{Colors.BOLD}Scaling Factors (from baseline):{Colors.RESET}")
            for snap in snapshots[1:]:
                x_factor = snap.domain_x / base.domain_x
                y_factor = snap.domain_y / base.domain_y
                z_factor = snap.domain_z / base.domain_z
                
                print(f"  • Node {snap.node_count}: "
                      f"X×{x_factor:.1f}, Y×{y_factor:.1f}, Z×{z_factor:.1f}")
        
        print()
    
    def print_step_log(self, step: int, node_count: int, active_dim: str, 
                       domain: Tuple[float, float, float],
                       cells: Tuple[int, int, int],
                       mpi: Tuple[int, int, int]):
        """
        Print detailed step-by-step log.
        
        Args:
            step: Step index
            node_count: Number of nodes
            active_dim: Active scaling dimension
            domain: Domain size tuple
            cells: Cell count tuple
            mpi: MPI decomposition tuple
        """
        if step == 0:
            print(f"\n{Colors.BOLD}{Colors.GREEN}Step {step}: Node {node_count} - BASELINE{Colors.RESET}")
            print(f"  Using exact values from run.yaml (no scaling)")
        else:
            print(f"\n{Colors.BOLD}Step {step}: Node {node_count} - Scale {Colors.RED}{active_dim}{Colors.RESET}")
        
        print(f"  Domain:  ({domain[0]:.1f}, {domain[1]:.1f}, {domain[2]:.1f})")
        print(f"  Cells:   ({cells[0]}, {cells[1]}, {cells[2]})")
        print(f"  MPI:     ({mpi[0]}, {mpi[1]}, {mpi[2]}) = {mpi[0]*mpi[1]*mpi[2]} procs")


def create_snapshot_from_job(job_config, active_dim: str = 'BASELINE') -> ScalingSnapshot:
    """
    Create ScalingSnapshot from JobConfig.
    
    Args:
        job_config: JobConfig object
        active_dim: Active dimension ('X', 'Y', 'Z', or 'BASELINE')
    
    Returns:
        ScalingSnapshot object
    """
    return ScalingSnapshot(
        node_count=job_config.num_nodes,
        num_procs=job_config.num_procs,
        active_dimension=active_dim,
        domain_x=job_config.domain_size[0] if job_config.domain_size else 0,
        domain_y=job_config.domain_size[1] if job_config.domain_size else 0,
        domain_z=job_config.domain_size[2] if job_config.domain_size else 0,
        cells_x=job_config.cell_count[0] if job_config.cell_count else 0,
        cells_y=job_config.cell_count[1] if job_config.cell_count else 0,
        cells_z=job_config.cell_count[2] if job_config.cell_count else 0,
        mpi_x=job_config.procs_decomposition[0],
        mpi_y=job_config.procs_decomposition[1],
        mpi_z=job_config.procs_decomposition[2],
        particles_x=1,  # Will be filled from config
        particles_y=1,
        particles_z=1
    )
