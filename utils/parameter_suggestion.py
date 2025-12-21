"""
Dynamic parameter suggestion for HPC simulations.

Automatically computes optimal grid points, particles per cell, and domain
decomposition based on detected hardware resources (cores, memory, GPUs).
"""

import logging
import math
from typing import Tuple, Dict, Any, Optional
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class SimulationParameters:
    """Suggested simulation parameters for a single node."""
    # Grid resolution
    nx: int
    ny: int
    nz: int
    
    # Particles per cell
    npcelx: int
    npcely: int
    npcelz: int
    
    # Domain decomposition
    core_x: int
    core_y: int
    core_z: int
    
    # Memory usage estimate
    estimated_memory_gb: float
    
    # Resource utilization
    cores_used: int
    cores_available: int
    memory_available_gb: float
    
    def __str__(self) -> str:
        return f"""
Suggested Simulation Parameters (Single Node):
  Grid Points:       nx={self.nx}, ny={self.ny}, nz={self.nz}
  Particles/Cell:    px={self.npcelx}, py={self.npcely}, pz={self.npcelz}
  Decomposition:     {self.core_x}×{self.core_y}×{self.core_z} = {self.cores_used} cores
  Memory Estimate:   {self.estimated_memory_gb:.1f} GB / {self.memory_available_gb:.1f} GB available
  Core Utilization:  {self.cores_used}/{self.cores_available} ({100*self.cores_used/self.cores_available:.1f}%)
"""


class ParameterSuggestion:
    """
    Dynamically suggest optimal simulation parameters based on hardware.
    
    This class computes grid resolution, particles per cell, and domain
    decomposition to maximize resource utilization on a single node.
    """
    
    def __init__(self, cores_per_node: int, memory_gb: float, 
                 gpus_per_node: int = 0, is_gpu_run: bool = False):
        """
        Initialize parameter suggestion engine.
        
        Args:
            cores_per_node: Number of CPU cores per node
            memory_gb: Available memory in GB
            gpus_per_node: Number of GPUs per node
            is_gpu_run: Whether this is a GPU-accelerated run
        """
        self.cores = cores_per_node
        self.memory_gb = memory_gb
        self.gpus = gpus_per_node
        self.is_gpu_run = is_gpu_run
        
        # Memory safety factor (use 80% of available memory)
        self.memory_safety_factor = 0.8
        
        # Typical memory per particle (bytes) - conservative estimate
        # Includes particle data, fields, and overhead
        self.bytes_per_particle = 200  # ~200 bytes per particle typical for PIC codes
        
        logger.debug(f"Parameter suggestion initialized: {cores_per_node} cores, "
                    f"{memory_gb:.1f} GB memory, {gpus_per_node} GPUs")
    
    def suggest_for_single_node(self, 
                                target_particles_per_cell: Optional[Tuple[int, int, int]] = None,
                                target_aspect_ratio: Tuple[float, float, float] = (1.0, 1.0, 1.0)
                                ) -> SimulationParameters:
        """
        Suggest optimal parameters for a single node.
        
        Args:
            target_particles_per_cell: Desired particles per cell (px, py, pz).
                                      If None, will compute optimal value.
            target_aspect_ratio: Desired domain aspect ratio (rx, ry, rz).
                                Default is cubic domain (1:1:1).
        
        Returns:
            SimulationParameters with suggested configuration
        """
        logger.info("Computing optimal parameters for single node...")
        
        # Step 1: Compute optimal domain decomposition
        core_x, core_y, core_z = self._compute_decomposition(
            self.cores, target_aspect_ratio
        )
        
        # Step 2: Determine particles per cell
        if target_particles_per_cell:
            npcelx, npcely, npcelz = target_particles_per_cell
        else:
            npcelx, npcely, npcelz = self._suggest_particles_per_cell()
        
        # Step 3: Compute grid resolution based on memory constraints
        nx, ny, nz = self._compute_grid_resolution(
            core_x, core_y, core_z,
            npcelx, npcely, npcelz,
            target_aspect_ratio
        )
        
        # Step 4: Estimate memory usage
        total_cells = nx * ny * nz
        total_particles = total_cells * npcelx * npcely * npcelz
        estimated_memory = (total_particles * self.bytes_per_particle) / (1024**3)
        
        params = SimulationParameters(
            nx=nx, ny=ny, nz=nz,
            npcelx=npcelx, npcely=npcely, npcelz=npcelz,
            core_x=core_x, core_y=core_y, core_z=core_z,
            estimated_memory_gb=estimated_memory,
            cores_used=core_x * core_y * core_z,
            cores_available=self.cores,
            memory_available_gb=self.memory_gb
        )
        
        logger.info(f"Suggested parameters: {nx}×{ny}×{nz} grid, "
                   f"{npcelx}×{npcely}×{npcelz} particles/cell, "
                   f"{core_x}×{core_y}×{core_z} decomposition")
        
        return params
    
    def _compute_decomposition(self, total_cores: int, 
                               aspect_ratio: Tuple[float, float, float]) -> Tuple[int, int, int]:
        """
        Compute optimal 3D domain decomposition.
        
        Args:
            total_cores: Total number of cores to decompose
            aspect_ratio: Target aspect ratio (rx, ry, rz)
        
        Returns:
            Tuple (core_x, core_y, core_z) such that product equals total_cores
        """
        rx, ry, rz = aspect_ratio
        
        # Find factorization closest to target aspect ratio
        best_decomp = None
        best_score = float('inf')
        
        for cx in range(1, total_cores + 1):
            if total_cores % cx != 0:
                continue
            
            remaining = total_cores // cx
            for cy in range(1, remaining + 1):
                if remaining % cy != 0:
                    continue
                
                cz = remaining // cy
                
                # Score based on deviation from aspect ratio
                current_ratio_x = cx / total_cores
                current_ratio_y = cy / total_cores
                current_ratio_z = cz / total_cores
                
                target_ratio_x = rx / (rx + ry + rz)
                target_ratio_y = ry / (rx + ry + rz)
                target_ratio_z = rz / (rx + ry + rz)
                
                score = (
                    abs(current_ratio_x - target_ratio_x) +
                    abs(current_ratio_y - target_ratio_y) +
                    abs(current_ratio_z - target_ratio_z)
                )
                
                if score < best_score:
                    best_score = score
                    best_decomp = (cx, cy, cz)
        
        if best_decomp is None:
            # Fallback: try to make it as cubic as possible
            cube_root = int(round(total_cores ** (1.0/3.0)))
            best_decomp = (cube_root, cube_root, total_cores // (cube_root * cube_root))
            logger.warning(f"Could not find perfect decomposition, using approximation: {best_decomp}")
        
        return best_decomp
    
    def _suggest_particles_per_cell(self) -> Tuple[int, int, int]:
        """
        Suggest number of particles per cell based on simulation type.
        
        For PIC codes, typical values are 8-64 particles per cell.
        We use conservative defaults that work well across applications.
        
        Returns:
            Tuple (npcelx, npcely, npcelz)
        """
        if self.is_gpu_run:
            # GPU runs can handle more particles efficiently
            return (4, 4, 4)  # 64 particles per cell
        else:
            # CPU runs use moderate particle count
            return (2, 2, 2)  # 8 particles per cell
    
    def _compute_grid_resolution(self, 
                                 core_x: int, core_y: int, core_z: int,
                                 npcelx: int, npcely: int, npcelz: int,
                                 aspect_ratio: Tuple[float, float, float]) -> Tuple[int, int, int]:
        """
        Compute grid resolution based on memory constraints.
        
        Args:
            core_x, core_y, core_z: Domain decomposition
            npcelx, npcely, npcelz: Particles per cell
            aspect_ratio: Target aspect ratio
        
        Returns:
            Tuple (nx, ny, nz) for grid resolution
        """
        # Available memory for particles (with safety factor)
        available_memory_bytes = self.memory_gb * self.memory_safety_factor * (1024**3)
        
        # Maximum total particles based on memory
        max_total_particles = int(available_memory_bytes / self.bytes_per_particle)
        
        # Particles per cell
        particles_per_cell = npcelx * npcely * npcelz
        
        # Maximum total cells
        max_total_cells = max_total_particles / particles_per_cell
        
        # Distribute cells according to aspect ratio
        rx, ry, rz = aspect_ratio
        volume_ratio = rx * ry * rz
        
        # Solve for cell counts maintaining aspect ratio
        # nx/ny/nz should be proportional to rx/ry/rz
        # nx * ny * nz = max_total_cells
        
        # For simplicity, solve assuming cubic root then scale
        base_cells = int((max_total_cells / volume_ratio) ** (1.0/3.0))
        
        nx = int(base_cells * rx)
        ny = int(base_cells * ry)
        nz = int(base_cells * rz)
        
        # Ensure cells are divisible by decomposition for even distribution
        nx = self._round_up_to_multiple(nx, core_x)
        ny = self._round_up_to_multiple(ny, core_y)
        nz = self._round_up_to_multiple(nz, core_z)
        
        # Ensure minimum grid size (at least 8 cells per dimension)
        nx = max(nx, core_x * 8)
        ny = max(ny, core_y * 8)
        nz = max(nz, core_z * 8)
        
        # Apply power-of-2 rounding for better performance (optional but recommended)
        nx = self._round_to_nice_number(nx)
        ny = self._round_to_nice_number(ny)
        nz = self._round_to_nice_number(nz)
        
        return (nx, ny, nz)
    
    def _round_up_to_multiple(self, value: int, multiple: int) -> int:
        """Round up to nearest multiple."""
        return ((value + multiple - 1) // multiple) * multiple
    
    def _round_to_nice_number(self, value: int) -> int:
        """
        Round to a 'nice' number (power of 2 or multiple of 64).
        
        This can improve cache performance and FFT efficiency.
        """
        # Find nearest power of 2
        power_of_2 = 2 ** int(math.log2(value))
        
        # Find nearest multiple of 64
        multiple_64 = ((value + 31) // 64) * 64
        
        # Choose whichever is closer but not smaller than original
        candidates = [power_of_2, power_of_2 * 2, multiple_64]
        valid = [c for c in candidates if c >= value]
        
        if valid:
            return min(valid, key=lambda x: abs(x - value))
        else:
            return value
    
    def suggest_for_weak_scaling(self, 
                                 num_nodes: int,
                                 base_params: Optional[SimulationParameters] = None
                                 ) -> SimulationParameters:
        """
        Suggest parameters for weak scaling test.
        
        For weak scaling, problem size increases proportionally with nodes,
        maintaining constant work per process.
        
        Args:
            num_nodes: Number of nodes for this scaling point
            base_params: Base parameters for single node (if None, will compute)
        
        Returns:
            SimulationParameters scaled for num_nodes
        """
        if base_params is None:
            base_params = self.suggest_for_single_node()
        
        # Scale factor for weak scaling (cubic root for 3D)
        scale_factor = num_nodes
        linear_scale = scale_factor ** (1.0 / 3.0)
        
        # Scale grid resolution
        nx = int(round(base_params.nx * linear_scale))
        ny = int(round(base_params.ny * linear_scale))
        nz = int(round(base_params.nz * linear_scale))
        
        # Particles per cell remain constant
        npcelx = base_params.npcelx
        npcely = base_params.npcely
        npcelz = base_params.npcelz
        
        # Scale decomposition
        total_cores = self.cores * num_nodes
        aspect_ratio = (
            base_params.core_x / (base_params.core_x + base_params.core_y + base_params.core_z),
            base_params.core_y / (base_params.core_x + base_params.core_y + base_params.core_z),
            base_params.core_z / (base_params.core_x + base_params.core_y + base_params.core_z)
        )
        
        core_x, core_y, core_z = self._compute_decomposition(total_cores, aspect_ratio)
        
        # Memory estimate
        total_cells = nx * ny * nz
        total_particles = total_cells * npcelx * npcely * npcelz
        estimated_memory = (total_particles * self.bytes_per_particle) / (1024**3)
        
        params = SimulationParameters(
            nx=nx, ny=ny, nz=nz,
            npcelx=npcelx, npcely=npcely, npcelz=npcelz,
            core_x=core_x, core_y=core_y, core_z=core_z,
            estimated_memory_gb=estimated_memory,
            cores_used=total_cores,
            cores_available=self.cores * num_nodes,
            memory_available_gb=self.memory_gb * num_nodes
        )
        
        return params
    
    def get_scaling_suggestions(self, 
                               max_nodes: int,
                               scaling_type: str = "weak") -> Dict[int, SimulationParameters]:
        """
        Get parameter suggestions for a complete scaling test.
        
        Args:
            max_nodes: Maximum number of nodes
            scaling_type: "weak" or "strong"
        
        Returns:
            Dictionary mapping node count to suggested parameters
        """
        suggestions = {}
        
        # Get base parameters for single node
        base_params = self.suggest_for_single_node()
        suggestions[1] = base_params
        
        # Generate power-of-2 sequence
        node_sequence = []
        n = 2
        while n <= max_nodes:
            node_sequence.append(n)
            n *= 2
        
        # Add max_nodes if not already included
        if max_nodes not in node_sequence and max_nodes > 1:
            node_sequence.append(max_nodes)
        
        for num_nodes in sorted(node_sequence):
            if scaling_type.lower() == "weak":
                suggestions[num_nodes] = self.suggest_for_weak_scaling(num_nodes, base_params)
            else:  # strong scaling
                # For strong scaling, grid stays constant, only decomposition changes
                total_cores = self.cores * num_nodes
                core_x, core_y, core_z = self._compute_decomposition(
                    total_cores, 
                    (base_params.core_x, base_params.core_y, base_params.core_z)
                )
                
                suggestions[num_nodes] = SimulationParameters(
                    nx=base_params.nx,
                    ny=base_params.ny,
                    nz=base_params.nz,
                    npcelx=base_params.npcelx,
                    npcely=base_params.npcely,
                    npcelz=base_params.npcelz,
                    core_x=core_x,
                    core_y=core_y,
                    core_z=core_z,
                    estimated_memory_gb=base_params.estimated_memory_gb,
                    cores_used=total_cores,
                    cores_available=self.cores * num_nodes,
                    memory_available_gb=self.memory_gb * num_nodes
                )
        
        return suggestions


def suggest_parameters(cores_per_node: int, 
                      memory_gb: float,
                      gpus_per_node: int = 0,
                      is_gpu_run: bool = False,
                      particles_per_cell: Optional[Tuple[int, int, int]] = None) -> SimulationParameters:
    """
    Convenience function to get parameter suggestions.
    
    Args:
        cores_per_node: CPU cores available per node
        memory_gb: Memory available in GB
        gpus_per_node: Number of GPUs per node
        is_gpu_run: Whether this is GPU-accelerated
        particles_per_cell: Optional user-specified particles per cell
    
    Returns:
        SimulationParameters with optimal configuration
    """
    suggester = ParameterSuggestion(cores_per_node, memory_gb, gpus_per_node, is_gpu_run)
    return suggester.suggest_for_single_node(target_particles_per_cell=particles_per_cell)
