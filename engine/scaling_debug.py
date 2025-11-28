"""
Debug-enabled scaling engine with comprehensive logging.

This module extends the base ScalingEngine with detailed step-by-step logging
for troubleshooting and validation of weak scaling configurations.

Usage:
    from engine.scaling_debug import DebugScalingEngine
    
    engine = DebugScalingEngine(config, resource_config, debug_level='VERBOSE')
    jobs = engine.generate_job_configs()
"""

import logging
from typing import List, Tuple
from enum import Enum

from engine.scaling import ScalingEngine, WeakScalingState
from core.config import ScalingConfig, JobConfig, ResourceConfig


class DebugLevel(str, Enum):
    """Debug verbosity levels."""
    MINIMAL = "minimal"      # Only key results
    STANDARD = "standard"    # Standard debug output
    VERBOSE = "verbose"      # Full step-by-step details
    TRACE = "trace"          # Every calculation logged


logger = logging.getLogger(__name__)


class DebugScalingEngine(ScalingEngine):
    """
    Scaling engine with comprehensive debug logging.
    
    Extends base ScalingEngine to provide detailed logging of:
    - Active dimension selection
    - Scale count calculations
    - Multiplier application
    - MPI correction steps
    - Domain/grid/decomposition evolution
    """
    
    def __init__(
        self,
        scaling_config: ScalingConfig,
        resource_config: ResourceConfig,
        debug_level: str = 'STANDARD'
    ):
        """
        Initialize debug scaling engine.
        
        Args:
            scaling_config: Scaling configuration from run.yaml
            resource_config: Resource configuration
            debug_level: Debug verbosity ('minimal', 'standard', 'verbose', 'trace')
        """
        super().__init__(scaling_config, resource_config, debug=True)
        
        self.debug_level = DebugLevel(debug_level.lower())
        
        # Configure logging based on debug level
        if self.debug_level == DebugLevel.TRACE:
            logger.setLevel(logging.DEBUG)
            logging.basicConfig(
                format='%(asctime)s [%(levelname)s] %(message)s',
                level=logging.DEBUG
            )
        else:
            logger.setLevel(logging.INFO)
        
        logger.info(f"╔{'═'*78}╗")
        logger.info(f"║ DEBUG SCALING ENGINE - {self.debug_level.value.upper()} MODE".ljust(78) + " ║")
        logger.info(f"╚{'═'*78}╝")
    
    def _apply_weak_scaling_step(
        self,
        step_index: int,
        num_nodes: int,
        base_domain: Tuple[float, float, float],
        base_cells: Tuple[int, int, int],
        base_procs: Tuple[int, int, int],
        base_particles: Tuple[int, int, int],
        factor: float,
        scaling_dims: int
    ) -> WeakScalingState:
        """
        Apply weak scaling step with detailed logging.
        
        Overrides parent method to add comprehensive debug output.
        """
        logger.info(f"\n{'#'*80}")
        logger.info(f"# STEP {step_index} → NODE {num_nodes}")
        logger.info(f"{'#'*80}")
        
        # Determine active dimension
        active_dim, scale_counts = self._get_active_dimension_and_counts(
            step_index, scaling_dims
        )
        x_count, y_count, z_count = scale_counts
        
        logger.info(f"\n┌─ DIMENSION SELECTION ─────────────────────────────────────")
        logger.info(f"│ Scaling mode: {scaling_dims}D")
        logger.info(f"│ Step index:   {step_index}")
        logger.info(f"│ Active dim:   {active_dim}")
        logger.info(f"└────────────────────────────────────────────────────────────")
        
        # Calculate multipliers
        x_mult = factor ** x_count
        y_mult = factor ** y_count
        z_mult = factor ** z_count
        
        logger.info(f"\n┌─ CUMULATIVE SCALE COUNTS ─────────────────────────────────")
        logger.info(f"│ X has scaled {x_count} times → multiplier = {factor}^{x_count} = {x_mult:.3f}")
        logger.info(f"│ Y has scaled {y_count} times → multiplier = {factor}^{y_count} = {y_mult:.3f}")
        logger.info(f"│ Z has scaled {z_count} times → multiplier = {factor}^{z_count} = {z_mult:.3f}")
        logger.info(f"└────────────────────────────────────────────────────────────")
        
        if self.debug_level in [DebugLevel.VERBOSE, DebugLevel.TRACE]:
            logger.info(f"\n┌─ BASELINE VALUES (from run.yaml) ─────────────────────────")
            logger.info(f"│ Domain: Lx={base_domain[0]}, Ly={base_domain[1]}, Lz={base_domain[2]}")
            logger.info(f"│ Cells:  nx={base_cells[0]}, ny={base_cells[1]}, nz={base_cells[2]}")
            logger.info(f"│ Procs:  px={base_procs[0]}, py={base_procs[1]}, pz={base_procs[2]}")
            logger.info(f"└────────────────────────────────────────────────────────────")
        
        # Apply scaling
        domain_scaled = (
            base_domain[0] * x_mult,
            base_domain[1] * y_mult,
            base_domain[2] * z_mult
        )
        
        cells_scaled = (
            int(round(base_cells[0] * x_mult)),
            int(round(base_cells[1] * y_mult)),
            int(round(base_cells[2] * z_mult))
        )
        
        procs_scaled = (
            int(base_procs[0] * x_mult),
            int(base_procs[1] * y_mult),
            int(base_procs[2] * z_mult)
        )
        
        logger.info(f"\n┌─ SCALED VALUES (before MPI correction) ───────────────────")
        logger.info(f"│ Domain: Lx={domain_scaled[0]:.2f} ({base_domain[0]:.2f}×{x_mult:.2f})")
        logger.info(f"│         Ly={domain_scaled[1]:.2f} ({base_domain[1]:.2f}×{y_mult:.2f})")
        logger.info(f"│         Lz={domain_scaled[2]:.2f} ({base_domain[2]:.2f}×{z_mult:.2f})")
        logger.info(f"│")
        logger.info(f"│ Cells:  nx={cells_scaled[0]} ({base_cells[0]}×{x_mult:.2f})")
        logger.info(f"│         ny={cells_scaled[1]} ({base_cells[1]}×{y_mult:.2f})")
        logger.info(f"│         nz={cells_scaled[2]} ({base_cells[2]}×{z_mult:.2f})")
        logger.info(f"│")
        logger.info(f"│ Procs:  px={procs_scaled[0]} ({base_procs[0]}×{x_mult:.2f})")
        logger.info(f"│         py={procs_scaled[1]} ({base_procs[1]}×{y_mult:.2f})")
        logger.info(f"│         pz={procs_scaled[2]} ({base_procs[2]}×{z_mult:.2f})")
        logger.info(f"└────────────────────────────────────────────────────────────")
        
        # MPI correction
        required_total_procs = num_nodes * self.resource_config.procs_per_node
        computed_total_procs = procs_scaled[0] * procs_scaled[1] * procs_scaled[2]
        
        logger.info(f"\n┌─ MPI RANK VALIDATION ─────────────────────────────────────")
        logger.info(f"│ Required: {num_nodes} nodes × {self.resource_config.procs_per_node} procs/node = {required_total_procs}")
        logger.info(f"│ Computed: {procs_scaled[0]}×{procs_scaled[1]}×{procs_scaled[2]} = {computed_total_procs}")
        
        if computed_total_procs != required_total_procs:
            logger.info(f"│ ⚠ MISMATCH DETECTED → Correcting {active_dim} dimension")
            procs_corrected = self._correct_mpi_decomposition(
                procs_scaled,
                active_dim,
                required_total_procs
            )
            logger.info(f"│ Correction: {procs_scaled} → {procs_corrected}")
            procs_scaled = procs_corrected
            computed_total_procs = procs_scaled[0] * procs_scaled[1] * procs_scaled[2]
            logger.info(f"│ ✓ Corrected: {procs_scaled[0]}×{procs_scaled[1]}×{procs_scaled[2]} = {computed_total_procs}")
        else:
            logger.info(f"│ ✓ No correction needed")
        
        logger.info(f"└────────────────────────────────────────────────────────────")
        
        # Create state
        state = WeakScalingState(
            node_count=num_nodes,
            step_index=step_index,
            active_dimension=active_dim,
            domain=domain_scaled,
            cells=cells_scaled,
            procs=procs_scaled,
            particles_per_cell=base_particles,
            total_procs=computed_total_procs
        )
        
        logger.info(f"\n┌─ FINAL CONFIGURATION ─────────────────────────────────────")
        logger.info(f"│ Node:       {num_nodes}")
        logger.info(f"│ Domain:     ({domain_scaled[0]:.2f}, {domain_scaled[1]:.2f}, {domain_scaled[2]:.2f})")
        logger.info(f"│ Cells:      ({cells_scaled[0]}, {cells_scaled[1]}, {cells_scaled[2]})")
        logger.info(f"│ MPI decomp: ({procs_scaled[0]}, {procs_scaled[1]}, {procs_scaled[2]}) = {computed_total_procs}")
        logger.info(f"│ Particles:  ({base_particles[0]}, {base_particles[1]}, {base_particles[2]}) [CONSTANT]")
        logger.info(f"└────────────────────────────────────────────────────────────")
        
        if self.debug_level == DebugLevel.TRACE:
            # Resolution check
            res_x = cells_scaled[0] / procs_scaled[0] if procs_scaled[0] > 0 else 0
            res_y = cells_scaled[1] / procs_scaled[1] if procs_scaled[1] > 0 else 0
            res_z = cells_scaled[2] / procs_scaled[2] if procs_scaled[2] > 0 else 0
            
            logger.info(f"\n┌─ RESOLUTION PER RANK ─────────────────────────────────────")
            logger.info(f"│ X: {cells_scaled[0]} / {procs_scaled[0]} = {res_x:.2f} cells/rank")
            logger.info(f"│ Y: {cells_scaled[1]} / {procs_scaled[1]} = {res_y:.2f} cells/rank")
            logger.info(f"│ Z: {cells_scaled[2]} / {procs_scaled[2]} = {res_z:.2f} cells/rank")
            logger.info(f"└────────────────────────────────────────────────────────────")
        
        return state
    
    def _correct_mpi_decomposition(
        self,
        procs: Tuple[int, int, int],
        active_dim: str,
        required_total: int
    ) -> Tuple[int, int, int]:
        """
        Correct MPI decomposition with detailed logging.
        
        Overrides parent method to add debug output.
        """
        px, py, pz = procs
        
        if self.debug_level == DebugLevel.TRACE:
            logger.debug(f"  ┌─ MPI CORRECTION DETAILS")
            logger.debug(f"  │ Input:    ({px}, {py}, {pz})")
            logger.debug(f"  │ Active:   {active_dim}")
            logger.debug(f"  │ Required: {required_total}")
        
        # Call parent implementation
        result = super()._correct_mpi_decomposition(procs, active_dim, required_total)
        
        if self.debug_level == DebugLevel.TRACE:
            logger.debug(f"  │ Output:   {result}")
            logger.debug(f"  └─")
        
        return result


def main():
    """Example usage of debug scaling engine."""
    from core.types import ScalingType
    
    # Example configuration
    config = ScalingConfig(
        scaling_type=ScalingType.WEAK,
        max_nodes=8,
        scaling_factor=2.0,
        scaling_dimensions=2,
        initial_procs=(14, 8, 1),
        initial_domain=(84.0, 48.0, 1.0),
        initial_cells=(840, 480, 1),
        particles_per_cell=(20, 20, 1)
    )
    
    resource_config = ResourceConfig(procs_per_node=112)
    
    # Create debug engine
    engine = DebugScalingEngine(config, resource_config, debug_level='VERBOSE')
    
    # Generate configurations
    jobs = engine.generate_job_configs()
    
    print(f"\n{'='*80}")
    print(f"GENERATED {len(jobs)} JOB CONFIGURATIONS")
    print(f"{'='*80}")
    for job in jobs:
        print(f"  {job.job_id}: {job.num_nodes} nodes, {job.num_procs} procs, decomp={job.procs_decomposition}")


if __name__ == '__main__':
    main()
