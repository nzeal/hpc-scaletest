"""
Directional Strong Scaling Algorithm

Implements intelligent directional doubling for strong scaling:
- 2D: Double in X → Y → X → Y → ... pattern
- 3D: Double in X → Y → Z → X → Y → Z → ... pattern

This ensures systematic exploration of the MPI topology space while
maintaining proper divisibility constraints.
"""

import logging
from typing import List, Tuple, Dict, Optional

logger = logging.getLogger(__name__)


class DirectionalStrongScaling:
    """
    Generate strong scaling sequence by directional doubling.
    
    For 2D scaling with base (14, 8, 1):
        Step 0: 14×8×1 = 112       (baseline)
        Step 1: 28×8×1 = 224       (double X: 14→28)
        Step 2: 28×16×1 = 448      (double Y: 8→16)
        Step 3: 56×16×1 = 896      (double X: 28→56)
        Step 4: 56×32×1 = 1792     (double Y: 16→32)
        Step 5: 112×32×1 = 3584    (double X: 56→112)
        Step 6: 112×64×1 = 7168    (double Y: 32→64)
        ...
    """
    
    def __init__(self, 
                 initial_procs: Tuple[int, int, int],
                 grid_dims: Tuple[int, int, int],
                 scaling_factor: float = 2.0,
                 scaling_dimensions: int = 2,
                 max_steps: int = 10):
        """
        Initialize directional scaling generator.
        
        Args:
            initial_procs: Base MPI decomposition (px0, py0, pz0)
            grid_dims: Fixed grid dimensions (nx, ny, nz)
            scaling_factor: Scaling factor (typically 2.0 for doubling)
            scaling_dimensions: Number of dimensions to scale (2 or 3)
            max_steps: Maximum number of scaling steps
        """
        self.px0, self.py0, self.pz0 = initial_procs
        self.nx, self.ny, self.nz = grid_dims
        self.scaling_factor = scaling_factor
        self.scaling_dims = scaling_dimensions
        self.max_steps = max_steps
        
        # Force pz=1 for 2D scaling
        if scaling_dimensions == 2:
            self.pz0 = 1
        
        logger.info(f"Initialized DirectionalStrongScaling:")
        logger.info(f"  Base decomposition: {self.px0}×{self.py0}×{self.pz0}")
        logger.info(f"  Grid (FIXED): {self.nx}×{self.ny}×{self.nz}")
        logger.info(f"  Scaling: {scaling_dimensions}D with factor {scaling_factor}")
    
    def generate_sequence(self) -> List[Dict]:
        """
        Generate complete scaling sequence with validation.
        
        Returns:
            list: Sequence of scaling configurations, each containing:
                - step: Step number
                - procs: (px, py, pz) tuple
                - total_ranks: Total MPI ranks
                - nodes_112: Nodes needed (112 cores/node)
                - valid: Whether decomposition divides grid evenly
                - cells_per_proc: Grid cells per process
                - direction: Which direction was doubled
        """
        sequence = []
        
        # Step 0: Baseline
        px, py, pz = self.px0, self.py0, self.pz0
        self._add_step(sequence, 0, px, py, pz, "baseline")
        
        # Generate scaling steps
        for step in range(1, self.max_steps + 1):
            # Determine which direction to scale
            if self.scaling_dims == 2:
                # 2D: Alternate X → Y → X → Y → ...
                direction = 'X' if (step % 2 == 1) else 'Y'
            else:
                # 3D: Cycle X → Y → Z → X → Y → Z → ...
                cycle = (step - 1) % 3
                direction = ['X', 'Y', 'Z'][cycle]
            
            # Apply scaling in the determined direction
            if direction == 'X':
                px = int(px * self.scaling_factor)
            elif direction == 'Y':
                py = int(py * self.scaling_factor)
            elif direction == 'Z':
                pz = int(pz * self.scaling_factor)
            
            self._add_step(sequence, step, px, py, pz, direction)
        
        return sequence
    
    def _add_step(self, sequence: List[Dict], step: int, 
                  px: int, py: int, pz: int, direction: str):
        """Add a scaling step to the sequence with validation."""
        total_ranks = px * py * pz
        nodes_112 = (total_ranks + 111) // 112  # Ceiling division
        
        # Check divisibility
        x_divisible = (self.nx % px == 0)
        y_divisible = (self.ny % py == 0)
        z_divisible = (self.nz % pz == 0) if pz > 1 else True
        
        valid = x_divisible and y_divisible and z_divisible
        
        # Calculate cells per process
        cells_per_proc = (self.nx * self.ny * self.nz) / total_ranks
        
        # Calculate actual grid cells per direction
        cells_x = self.nx / px if px > 0 else 0
        cells_y = self.ny / py if py > 0 else 0
        cells_z = self.nz / pz if pz > 0 else 0
        
        config = {
            'step': step,
            'procs': (px, py, pz),
            'total_ranks': total_ranks,
            'nodes_112': nodes_112,
            'valid': valid,
            'cells_per_proc': cells_per_proc,
            'cells_per_dir': (cells_x, cells_y, cells_z),
            'direction': direction,
            'divisibility': {
                'x': x_divisible,
                'y': y_divisible,
                'z': z_divisible
            }
        }
        
        sequence.append(config)
        
        # Log the step
        status = "✓" if valid else "❌"
        logger.info(f"{status} Step {step}: {px}×{py}×{pz} = {total_ranks} ranks "
                   f"(direction: {direction}, nodes: {nodes_112})")
        
        if not valid:
            logger.warning(f"  Invalid decomposition: nx/px={self.nx}/{px}={cells_x:.2f}, "
                          f"ny/py={self.ny}/{py}={cells_y:.2f}, nz/pz={self.nz}/{pz}={cells_z:.2f}")
    
    def print_sequence(self, sequence: List[Dict]):
        """Print formatted scaling sequence."""
        print()
        print("="*85)
        print(f"Directional Strong Scaling Sequence ({self.scaling_dims}D)")
        print("="*85)
        print()
        print(f"Base configuration: {self.px0}×{self.py0}×{self.pz0} = {self.px0*self.py0*self.pz0} ranks")
        print(f"Grid (FIXED): {self.nx}×{self.ny}×{self.nz} = {self.nx*self.ny*self.nz:,} cells")
        print(f"Scaling: {self.scaling_factor}x in each direction")
        print()
        print(f"{'Step':<6} {'Direction':<10} {'Decomposition':<15} {'Ranks':<8} "
              f"{'Nodes':<7} {'Valid':<7} {'Cells/Proc':<12}")
        print("-"*85)
        
        for config in sequence:
            step = config['step']
            direction = config['direction']
            px, py, pz = config['procs']
            decomp = f"{px}×{py}×{pz}"
            ranks = config['total_ranks']
            nodes = config['nodes_112']
            valid = "✓" if config['valid'] else "❌"
            cells = config['cells_per_proc']
            
            print(f"{step:<6} {direction:<10} {decomp:<15} {ranks:<8} {nodes:<7} "
                  f"{valid:<7} {cells:<12.1f}")
        
        print("-"*85)
        print()
        
        # Summary
        valid_count = sum(1 for c in sequence if c['valid'])
        invalid_count = len(sequence) - valid_count
        
        print(f"Summary:")
        print(f"  Total steps: {len(sequence)}")
        print(f"  Valid configurations: {valid_count}")
        print(f"  Invalid configurations: {invalid_count}")
        
        if invalid_count > 0:
            print()
            print(f"❌ Invalid steps:")
            for config in sequence:
                if not config['valid']:
                    step = config['step']
                    px, py, pz = config['procs']
                    ranks = config['total_ranks']
                    nodes = config['nodes_112']
                    cells_x, cells_y, cells_z = config['cells_per_dir']
                    
                    print(f"  Step {step}: {px}×{py}×{pz} = {ranks} ranks ({nodes} nodes)")
                    print(f"    Problem: ", end="")
                    issues = []
                    if not config['divisibility']['x']:
                        issues.append(f"nx/px = {self.nx}/{px} = {cells_x:.2f}")
                    if not config['divisibility']['y']:
                        issues.append(f"ny/py = {self.ny}/{py} = {cells_y:.2f}")
                    if not config['divisibility']['z'] and pz > 1:
                        issues.append(f"nz/pz = {self.nz}/{pz} = {cells_z:.2f}")
                    print(", ".join(issues))
        
        print()
        print("="*85)
        print()
    
    def get_valid_sequence(self, sequence: List[Dict]) -> List[Dict]:
        """Filter to only valid configurations."""
        return [c for c in sequence if c['valid']]
    
    def suggest_alternatives(self, invalid_config: Dict) -> List[Tuple]:
        """Suggest alternative decompositions near an invalid configuration."""
        target_ranks = invalid_config['total_ranks']
        
        # Search for nearby rank counts with valid decompositions
        alternatives = []
        for ranks in range(int(target_ranks * 0.85), int(target_ranks * 1.15)):
            # Try factoring this rank count
            decomp = self._find_valid_decomposition(ranks)
            if decomp:
                px, py, pz = decomp
                nodes_112 = (ranks + 111) // 112
                diff = ranks - target_ranks
                diff_pct = (diff / target_ranks) * 100
                alternatives.append((ranks, decomp, nodes_112, diff, diff_pct))
        
        # Sort by proximity
        alternatives.sort(key=lambda x: abs(x[3]))
        return alternatives[:5]
    
    def _find_valid_decomposition(self, total_ranks: int) -> Optional[Tuple[int, int, int]]:
        """Find valid decomposition for given rank count."""
        # Get divisors
        rank_divisors = self._get_divisors(total_ranks)
        nx_divisors = set(self._get_divisors(self.nx))
        ny_divisors = set(self._get_divisors(self.ny))
        
        best = None
        best_score = float('inf')
        
        for px in rank_divisors:
            if self.scaling_dims == 2:
                py = total_ranks // px
                pz = 1
                
                if px * py != total_ranks:
                    continue
                
                if px in nx_divisors and py in ny_divisors:
                    # Calculate aspect ratio score
                    grid_aspect = self.nx / self.ny if self.ny > 0 else float('inf')
                    decomp_aspect = px / py if py > 0 else float('inf')
                    score = abs(grid_aspect - decomp_aspect)
                    
                    if score < best_score:
                        best_score = score
                        best = (px, py, pz)
        
        return best
    
    def _get_divisors(self, n: int) -> List[int]:
        """Get all divisors of n."""
        if n <= 0:
            return [1]
        divisors = []
        for i in range(1, int(n**0.5) + 1):
            if n % i == 0:
                divisors.append(i)
                if i != n // i:
                    divisors.append(n // i)
        return sorted(divisors)


if __name__ == "__main__":
    # Test directional scaling
    logging.basicConfig(level=logging.INFO)
    
    print("Testing Directional Strong Scaling Algorithm")
    print()
    
    # Test case: User's configuration (840×480, base 14×8×1)
    print("Test: User's Configuration (840×480 grid, 14×8×1 base)")
    scaler = DirectionalStrongScaling(
        initial_procs=(14, 8, 1),
        grid_dims=(840, 480, 1),
        scaling_factor=2.0,
        scaling_dimensions=2,
        max_steps=10
    )
    
    sequence = scaler.generate_sequence()
    scaler.print_sequence(sequence)
    
    # Show valid-only sequence
    valid_seq = scaler.get_valid_sequence(sequence)
    print("Valid configurations only:")
    for config in valid_seq:
        px, py, pz = config['procs']
        ranks = config['total_ranks']
        nodes = config['nodes_112']
        print(f"  Step {config['step']}: {px}×{py}×{pz} = {ranks} ranks ({nodes} nodes)")
