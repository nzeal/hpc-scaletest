"""
Smart Directional Strong Scaling with Hardware Constraints

Handles cases where perfect directional doubling is mathematically impossible.
Provides intelligent fallback strategies while maintaining strong scaling principles.
"""

import logging
from typing import List, Tuple, Dict, Optional
from math import gcd

logger = logging.getLogger(__name__)


class SmartDirectionalScaling:
    """
    Intelligent directional strong scaling that adapts to hardware constraints.
    
    Key features:
    1. Attempts directional doubling (X→Y or X→Y→Z pattern)
    2. Falls back to nearest valid decomposition if doubling fails
    3. Respects fixed grid dimensions and hardware cores/node
    4. Provides clear explanations for all decisions
    """
    
    def __init__(self,
                 initial_procs: Tuple[int, int, int],
                 grid_dims: Tuple[int, int, int],
                 procs_per_node: int,
                 max_nodes: int,  # REQUIRED - must be provided from YAML
                 scaling_factor: float = 2.0,
                 scaling_dimensions: int = 2):
        """
        Initialize smart directional scaling.
        
        Args:
            initial_procs: Base MPI decomposition (px0, py0, pz0)
            grid_dims: Fixed grid dimensions (nx, ny, nz)
            procs_per_node: Hardware constraint (from YAML or detection)
            max_nodes: Maximum number of nodes (REQUIRED - from YAML)
            scaling_factor: Scaling factor (typically 2.0)
            scaling_dimensions: Dimensions to scale (2 or 3)
        """
        self.px0, self.py0, self.pz0 = initial_procs
        self.nx, self.ny, self.nz = grid_dims
        self.procs_per_node = procs_per_node
        self.scaling_factor = scaling_factor
        self.scaling_dims = scaling_dimensions
        self.max_nodes = max_nodes
        
        # Force pz=1 for 2D
        if scaling_dimensions == 2:
            self.pz0 = 1
        
        # Precompute divisors for efficiency
        self.nx_divisors = set(self._get_divisors(self.nx))
        self.ny_divisors = set(self._get_divisors(self.ny))
        self.nz_divisors = set(self._get_divisors(self.nz))
        
        logger.info(f"Initialized SmartDirectionalScaling:")
        logger.info(f"  Base: {self.px0}×{self.py0}×{self.pz0} = {self.px0*self.py0*self.pz0}")
        logger.info(f"  Grid: {self.nx}×{self.ny}×{self.nz}")
        logger.info(f"  Hardware: {procs_per_node} cores/node")
        logger.info(f"  Mode: {scaling_dimensions}D scaling")
    
    def generate_sequence(self) -> List[Dict]:
        """
        Generate smart scaling sequence.
        
        Returns:
            list: Scaling configurations with validation and fallback info
        """
        sequence = []
        current_px, current_py, current_pz = self.px0, self.py0, self.pz0
        
        # Step 0: Baseline
        self._add_configuration(sequence, 0, current_px, current_py, current_pz,
                               "baseline", "exact")
        
        step = 0
        while True:
            step += 1
            
            # Calculate target node count (power-of-2 sequence)
            target_nodes = 2 ** step
            if target_nodes > self.max_nodes:
                break
            
            target_ranks = target_nodes * self.procs_per_node
            
            # Determine next direction to scale
            if self.scaling_dims == 2:
                next_dir = 'X' if (step % 2 == 1) else 'Y'
            else:
                next_dir = ['X', 'Y', 'Z'][(step - 1) % 3]
            
            # Try directional doubling
            if next_dir == 'X':
                trial_px = int(current_px * self.scaling_factor)
                trial_py, trial_pz = current_py, current_pz
            elif next_dir == 'Y':
                trial_px = current_px
                trial_py = int(current_py * self.scaling_factor)
                trial_pz = current_pz
            else:  # Z
                trial_px, trial_py = current_px, current_py
                trial_pz = int(current_pz * self.scaling_factor)
            
            # Check if trial decomposition is valid
            if self._is_valid_decomposition(trial_px, trial_py, trial_pz):
                # Success - use directional doubling
                current_px, current_py, current_pz = trial_px, trial_py, trial_pz
                self._add_configuration(sequence, step, current_px, current_py, current_pz,
                                       next_dir, "exact")
            else:
                # Directional doubling fails - find best alternative
                alt_decomp = self._find_best_decomposition(target_ranks)
                
                if alt_decomp:
                    current_px, current_py, current_pz = alt_decomp
                    self._add_configuration(sequence, step, current_px, current_py, current_pz,
                                           next_dir, "fallback")
                else:
                    # No valid decomposition found - skip this node count
                    self._add_configuration(sequence, step, None, None, None,
                                           next_dir, "impossible", target_ranks)
        
        return sequence
    
    def _is_valid_decomposition(self, px: int, py: int, pz: int) -> bool:
        """Check if decomposition divides grid evenly."""
        return (px in self.nx_divisors and 
                py in self.ny_divisors and
                (pz in self.nz_divisors or pz == 1))
    
    def _find_best_decomposition(self, target_ranks: int) -> Optional[Tuple[int, int, int]]:
        """
        Find best valid decomposition for target rank count.
        
        Priority:
        1. Exact match with perfect divisibility
        2. Close match (within 10%) with perfect divisibility
        3. None if no good match exists
        """
        # Try exact match first
        exact = self._find_decomposition_for_ranks(target_ranks)
        if exact:
            return exact
        
        # Search nearby rank counts (±10%)
        min_ranks = int(target_ranks * 0.90)
        max_ranks = int(target_ranks * 1.10)
        
        best = None
        best_diff = float('inf')
        
        for ranks in range(min_ranks, max_ranks + 1):
            decomp = self._find_decomposition_for_ranks(ranks)
            if decomp:
                diff = abs(ranks - target_ranks)
                if diff < best_diff:
                    best_diff = diff
                    best = decomp
        
        return best
    
    def _find_decomposition_for_ranks(self, total_ranks: int) -> Optional[Tuple[int, int, int]]:
        """Find valid decomposition for exact rank count."""
        if self.scaling_dims == 2:
            return self._find_2d_decomposition(total_ranks)
        else:
            return self._find_3d_decomposition(total_ranks)
    
    def _find_2d_decomposition(self, ranks: int) -> Optional[Tuple[int, int]]:
        """Find px, py such that px×py=ranks with divisibility."""
        rank_divisors = self._get_divisors(ranks)
        
        best = None
        best_score = float('inf')
        
        for px in rank_divisors:
            py = ranks // px
            
            if px * py != ranks:
                continue
            
            if px in self.nx_divisors and py in self.ny_divisors:
                # Calculate aspect ratio score
                grid_aspect = self.nx / self.ny if self.ny > 0 else float('inf')
                decomp_aspect = px / py if py > 0 else float('inf')
                score = abs(grid_aspect - decomp_aspect)
                
                if score < best_score:
                    best_score = score
                    best = (px, py, 1)
        
        return best
    
    def _find_3d_decomposition(self, ranks: int) -> Optional[Tuple[int, int, int]]:
        """Find px, py, pz such that px×py×pz=ranks with divisibility."""
        # TODO: Implement full 3D search
        # For now, use 2D with pz=1
        decomp_2d = self._find_2d_decomposition(ranks)
        return decomp_2d if decomp_2d else None
    
    def _add_configuration(self, sequence: List[Dict], step: int,
                          px: Optional[int], py: Optional[int], pz: Optional[int],
                          direction: str, strategy: str, target_ranks: int = None):
        """Add configuration to sequence."""
        if px is None:
            # Impossible configuration
            config = {
                'step': step,
                'target_nodes': 2 ** step,
                'target_ranks': target_ranks,
                'procs': None,
                'valid': False,
                'strategy': strategy,
                'direction': direction,
                'reason': 'No valid decomposition exists'
            }
        else:
            total_ranks = px * py * pz
            nodes_needed = (total_ranks + self.procs_per_node - 1) // self.procs_per_node
            
            config = {
                'step': step,
                'target_nodes': 2 ** step,
                'actual_nodes': nodes_needed,
                'procs': (px, py, pz),
                'total_ranks': total_ranks,
                'valid': True,
                'strategy': strategy,
                'direction': direction,
                'cells_per_proc': (self.nx * self.ny * self.nz) / total_ranks,
                'cells_per_dir': (self.nx / px, self.ny / py, self.nz / pz)
            }
        
        sequence.append(config)
        
        # Logging
        if config['valid']:
            status = "✓ exact" if strategy == "exact" else "⚠ fallback"
            logger.info(f"{status} Step {step}: {px}×{py}×{pz} = {total_ranks} ranks "
                       f"({nodes_needed} nodes, dir: {direction})")
        else:
            logger.warning(f"❌ Step {step}: Target {target_ranks} ranks - IMPOSSIBLE")
    
    def print_sequence(self, sequence: List[Dict]):
        """Print formatted sequence with recommendations."""
        print()
        print("="*100)
        print(f"Smart Directional Strong Scaling ({self.scaling_dims}D)")
        print("="*100)
        print()
        print(f"Configuration:")
        print(f"  Grid (FIXED): {self.nx}×{self.ny}×{self.nz} = {self.nx*self.ny*self.nz:,} cells")
        print(f"  Base decomposition: {self.px0}×{self.py0}×{self.pz0} = {self.px0*self.py0*self.pz0} ranks")
        print(f"  Hardware: {self.procs_per_node} cores/node")
        print(f"  Scaling: {self.scaling_factor}x per step in {['X', 'Y'] if self.scaling_dims == 2 else ['X', 'Y', 'Z']} directions")
        print()
        
        print(f"{'Step':<6} {'Target':<8} {'Actual':<8} {'Decomposition':<16} {'Ranks':<8} {'Status':<12} {'Strategy':<10}")
        print("-"*100)
        
        valid_configs = []
        invalid_configs = []
        
        for config in sequence:
            step = config['step']
            target = config['target_nodes']
            
            if config['valid']:
                actual = config['actual_nodes']
                px, py, pz = config['procs']
                decomp = f"{px}×{py}×{pz}"
                ranks = config['total_ranks']
                status = "✓ Valid" if config['strategy'] == "exact" else "⚠ Adjusted"
                strategy = config['strategy'].capitalize()
                
                print(f"{step:<6} {target:<8} {actual:<8} {decomp:<16} {ranks:<8} {status:<12} {strategy:<10}")
                valid_configs.append(config)
            else:
                actual = "—"
                decomp = "—"
                ranks = config.get('target_ranks', '?')
                status = "❌ Skip"
                strategy = "None"
                
                print(f"{step:<6} {target:<8} {actual:<8} {decomp:<16} {ranks:<8} {status:<12} {strategy:<10}")
                invalid_configs.append(config)
        
        print("-"*100)
        print()
        
        # Recommendations
        print("="*100)
        print("RECOMMENDATIONS")
        print("="*100)
        print()
        
        if invalid_configs:
            print(f"⚠ {len(invalid_configs)} node count(s) have NO valid decomposition")
            print()
            print("For your YAML configuration, use:")
            print(f"  node_counts: {[c['actual_nodes'] for c in valid_configs]}")
            print()
        else:
            print("✅ All power-of-2 node counts have valid decompositions!")
            print()
        
        print("="*100)
        print()
    
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
    logging.basicConfig(level=logging.INFO)
    
    print("Testing Smart Directional Strong Scaling")
    print()
    
    # User's configuration
    scaler = SmartDirectionalScaling(
        initial_procs=(14, 8, 1),
        grid_dims=(840, 480, 1),
        procs_per_node=112,
        scaling_factor=2.0,
        scaling_dimensions=2,
        max_nodes=128
    )
    
    sequence = scaler.generate_sequence()
    scaler.print_sequence(sequence)
