"""
Node-Sequence-Based Strong Scaling (CORRECT Implementation)

Mathematical Logic:
1. Generate node sequence: nodes = [scaling_factor^i for i in range(max_steps)]
   - scaling_factor=2.0 → [1, 2, 4, 8, 16, 32, 64, 128]
   - scaling_factor=4.0 → [1, 4, 16, 64, 256]

2. For each node count:
   - Step 0: Use initial_procs EXACTLY as given (baseline)
   - Step 1+: Calculate target_ranks = nodes × procs_per_node
             Find best valid decomposition for target_ranks
             Prefer directional doubling pattern (X→Y→Z)

Key Insight: The NODE SEQUENCE is determined by scaling_factor,
             NOT by doubling the decomposition itself!
"""

import logging
from typing import List, Tuple, Dict, Optional
from math import ceil

logger = logging.getLogger(__name__)


class NodeSequenceStrongScaling:
    """
    Generate strong scaling based on node sequence with directional preference.
    
    This is the CORRECT approach:
    - Node sequence determined by: nodes = [scaling_factor^i]
    - For each node count, find valid decomposition
    - Prefer decompositions that follow X→Y→Z doubling pattern
    """
    
    def __init__(self,
                 initial_procs: Tuple[int, int, int],
                 grid_dims: Tuple[int, int, int],
                 procs_per_node: int,
                 scaling_factor: float = 2.0,
                 scaling_dimensions: int = 2,
                 max_nodes: int = 128):
        """
        Initialize node-sequence-based scaling.
        
        Args:
            initial_procs: Base MPI decomposition (px0, py0, pz0) - used EXACTLY for step 0
            grid_dims: Fixed grid dimensions (nx, ny, nz)
            procs_per_node: Hardware cores per node
            scaling_factor: Scaling factor (2.0 for doubling, 4.0 for quadrupling)
            scaling_dimensions: Dimensions to scale (2 or 3)
            max_nodes: Maximum nodes
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
        
        # Precompute divisors
        self.nx_divisors = set(self._get_divisors(self.nx))
        self.ny_divisors = set(self._get_divisors(self.ny))
        self.nz_divisors = set(self._get_divisors(self.nz))
        
        logger.info(f"Initialized NodeSequenceStrongScaling:")
        logger.info(f"  Base: {self.px0}×{self.py0}×{self.pz0} = {self.px0*self.py0*self.pz0}")
        logger.info(f"  Grid: {self.nx}×{self.ny}×{self.nz}")
        logger.info(f"  Hardware: {procs_per_node} cores/node")
        logger.info(f"  Scaling factor: {scaling_factor}")
    
    def generate_node_sequence(self) -> List[int]:
        """
        Generate node sequence based on scaling_factor.
        
        Returns:
            list: Node counts [scaling_factor^0, scaling_factor^1, ...]
        """
        sequence = []
        step = 0
        
        while True:
            nodes = int(self.scaling_factor ** step)
            if nodes > self.max_nodes:
                break
            sequence.append(nodes)
            step += 1
        
        logger.info(f"Generated node sequence: {sequence}")
        return sequence
    
    def generate_sequence(self) -> List[Dict]:
        """
        Generate complete scaling sequence.
        
        Returns:
            list: Scaling configurations
        """
        node_sequence = self.generate_node_sequence()
        configurations = []
        
        for step, nodes in enumerate(node_sequence):
            if step == 0:
                # Step 0: Use initial_procs EXACTLY
                px, py, pz = self.px0, self.py0, self.pz0
                total_ranks = px * py * pz
                
                # Verify it matches hardware
                expected_ranks = nodes * self.procs_per_node
                if total_ranks != expected_ranks:
                    logger.warning(f"Step 0: initial_procs ({total_ranks} ranks) "
                                 f"doesn't match {nodes} node(s) × {self.procs_per_node} cores/node "
                                 f"= {expected_ranks} ranks")
                
                config = {
                    'step': step,
                    'target_nodes': nodes,
                    'actual_nodes': nodes,
                    'procs': (px, py, pz),
                    'total_ranks': total_ranks,
                    'valid': self._is_valid_decomposition(px, py, pz),
                    'strategy': 'baseline',
                    'cells_per_proc': (self.nx * self.ny * self.nz) / total_ranks,
                }
                
                configurations.append(config)
                logger.info(f"✓ Step {step}: {px}×{py}×{pz} = {total_ranks} ranks "
                          f"({nodes} nodes) [baseline - user defined]")
            else:
                # Step 1+: Find best decomposition for target ranks
                target_ranks = nodes * self.procs_per_node
                
                # Find decomposition with directional preference
                decomp = self._find_best_decomposition_with_preference(
                    target_ranks, 
                    step,
                    configurations[-1]['procs']  # Previous decomposition for guidance
                )
                
                if decomp:
                    px, py, pz = decomp
                    actual_ranks = px * py * pz
                    actual_nodes = ceil(actual_ranks / self.procs_per_node)
                    
                    config = {
                        'step': step,
                        'target_nodes': nodes,
                        'actual_nodes': actual_nodes,
                        'procs': (px, py, pz),
                        'total_ranks': actual_ranks,
                        'valid': True,
                        'strategy': 'exact' if actual_ranks == target_ranks else 'approximate',
                        'cells_per_proc': (self.nx * self.ny * self.nz) / actual_ranks,
                        'deviation': ((actual_ranks - target_ranks) / target_ranks * 100) if target_ranks > 0 else 0
                    }
                    
                    status = "✓" if config['strategy'] == 'exact' else "⚠"
                    logger.info(f"{status} Step {step}: {px}×{py}×{pz} = {actual_ranks} ranks "
                              f"({actual_nodes} nodes, target: {nodes})")
                else:
                    # No valid decomposition found
                    config = {
                        'step': step,
                        'target_nodes': nodes,
                        'actual_nodes': None,
                        'procs': None,
                        'total_ranks': target_ranks,
                        'valid': False,
                        'strategy': 'impossible',
                        'reason': f'No valid decomposition for {target_ranks} ranks'
                    }
                    
                    logger.warning(f"❌ Step {step}: Target {nodes} nodes ({target_ranks} ranks) - "
                                 f"NO VALID DECOMPOSITION")
                
                configurations.append(config)
        
        return configurations
    
    def _find_best_decomposition_with_preference(self, 
                                                 target_ranks: int,
                                                 step: int,
                                                 prev_decomp: Tuple[int, int, int]) -> Optional[Tuple[int, int, int]]:
        """
        Find best decomposition for target ranks with directional preference.
        
        Strategy:
        1. Try exact match with directional doubling from previous step
        2. Try exact match with any valid decomposition
        3. Try approximate match (±10%) with directional preference
        4. Try approximate match with any valid decomposition
        """
        # Strategy 1: Exact match with directional doubling
        if prev_decomp:
            prev_px, prev_py, prev_pz = prev_decomp
            
            # Determine which direction to double based on step and scaling_dims
            if self.scaling_dims == 2:
                next_dir = 'X' if (step % 2 == 1) else 'Y'
            else:
                next_dir = ['X', 'Y', 'Z'][(step - 1) % 3]
            
            # Try directional doubling
            if next_dir == 'X':
                trial = (int(prev_px * self.scaling_factor), prev_py, prev_pz)
            elif next_dir == 'Y':
                trial = (prev_px, int(prev_py * self.scaling_factor), prev_pz)
            else:  # Z
                trial = (prev_px, prev_py, int(prev_pz * self.scaling_factor))
            
            px, py, pz = trial
            if px * py * pz == target_ranks and self._is_valid_decomposition(px, py, pz):
                logger.debug(f"  Found via directional doubling: {px}×{py}×{pz}")
                return trial
        
        # Strategy 2: Exact match with any decomposition
        exact = self._find_decomposition_for_ranks(target_ranks)
        if exact:
            logger.debug(f"  Found exact match: {exact}")
            return exact
        
        # Strategy 3: Approximate match (±10%)
        min_ranks = int(target_ranks * 0.90)
        max_ranks = int(target_ranks * 1.10)
        
        best = None
        best_score = float('inf')
        
        for ranks in range(min_ranks, max_ranks + 1):
            decomp = self._find_decomposition_for_ranks(ranks)
            if decomp:
                # Score based on proximity to target
                diff = abs(ranks - target_ranks)
                score = diff / target_ranks
                
                # Bonus for maintaining directional pattern
                if prev_decomp:
                    px, py, pz = decomp
                    prev_px, prev_py, prev_pz = prev_decomp
                    
                    # Check if this follows X→Y pattern
                    if self.scaling_dims == 2:
                        if step % 2 == 1:  # Should double X
                            if px > prev_px and py == prev_py:
                                score *= 0.8  # Prefer this
                        else:  # Should double Y
                            if py > prev_py and px == prev_px:
                                score *= 0.8
                
                if score < best_score:
                    best_score = score
                    best = decomp
        
        if best:
            logger.debug(f"  Found approximate match: {best} (deviation: {best_score*100:.1f}%)")
        
        return best
    
    def _is_valid_decomposition(self, px: int, py: int, pz: int) -> bool:
        """Check if decomposition divides grid evenly."""
        return (px in self.nx_divisors and 
                py in self.ny_divisors and
                (pz in self.nz_divisors or pz == 1))
    
    def _find_decomposition_for_ranks(self, total_ranks: int) -> Optional[Tuple[int, int, int]]:
        """Find valid decomposition for exact rank count."""
        if self.scaling_dims == 2:
            return self._find_2d_decomposition(total_ranks)
        else:
            return self._find_3d_decomposition(total_ranks)
    
    def _find_2d_decomposition(self, ranks: int) -> Optional[Tuple[int, int]]:
        """Find px, py such that px×py=ranks with perfect divisibility."""
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
        # For now, use 2D with pz=1
        decomp_2d = self._find_2d_decomposition(ranks)
        return decomp_2d if decomp_2d else None
    
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
    
    def print_sequence(self, configurations: List[Dict]):
        """Print formatted sequence."""
        print()
        print("="*100)
        print(f"Node-Sequence-Based Strong Scaling ({self.scaling_dims}D)")
        print("="*100)
        print()
        print(f"Configuration:")
        print(f"  Grid (FIXED): {self.nx}×{self.ny}×{self.nz} = {self.nx*self.ny*self.nz:,} cells")
        print(f"  Base: {self.px0}×{self.py0}×{self.pz0} = {self.px0*self.py0*self.pz0} ranks")
        print(f"  Hardware: {self.procs_per_node} cores/node")
        print(f"  Scaling factor: {self.scaling_factor} (nodes = {self.scaling_factor}^step)")
        print()
        
        # Show node sequence
        node_seq = [c['target_nodes'] for c in configurations if c.get('target_nodes')]
        print(f"Node sequence: {node_seq}")
        print()
        
        print(f"{'Step':<6} {'Target':<8} {'Actual':<8} {'Decomposition':<16} "
              f"{'Ranks':<8} {'Deviation':<12} {'Status':<10}")
        print("-"*100)
        
        for config in configurations:
            step = config['step']
            target = config['target_nodes']
            
            if config['valid']:
                actual = config['actual_nodes']
                px, py, pz = config['procs']
                decomp = f"{px}×{py}×{pz}"
                ranks = config['total_ranks']
                
                if 'deviation' in config:
                    dev = f"{config['deviation']:+.1f}%"
                else:
                    dev = "—"
                
                status = "✓ Valid" if config['strategy'] in ['baseline', 'exact'] else "⚠ Approx"
                
                print(f"{step:<6} {target:<8} {actual:<8} {decomp:<16} {ranks:<8} {dev:<12} {status:<10}")
            else:
                actual = "—"
                decomp = "—"
                ranks = config['total_ranks']
                dev = "—"
                status = "❌ Skip"
                
                print(f"{step:<6} {target:<8} {actual:<8} {decomp:<16} {ranks:<8} {dev:<12} {status:<10}")
        
        print("-"*100)
        print()
        
        # Summary
        valid_count = sum(1 for c in configurations if c['valid'])
        print(f"Summary: {valid_count}/{len(configurations)} valid configurations")
        print()
        print("="*100)
        print()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    print("Testing Node-Sequence-Based Strong Scaling")
    print()
    
    # Test 1: scaling_factor = 2.0
    print("="*80)
    print("Test 1: scaling_factor = 2.0 (doubling)")
    print("="*80)
    scaler1 = NodeSequenceStrongScaling(
        initial_procs=(14, 8, 1),
        grid_dims=(840, 480, 1),
        procs_per_node=112,
        scaling_factor=2.0,
        scaling_dimensions=2,
        max_nodes=128
    )
    seq1 = scaler1.generate_sequence()
    scaler1.print_sequence(seq1)
    
    # Test 2: scaling_factor = 4.0
    print()
    print("="*80)
    print("Test 2: scaling_factor = 4.0 (quadrupling)")
    print("="*80)
    scaler2 = NodeSequenceStrongScaling(
        initial_procs=(14, 8, 1),
        grid_dims=(840, 480, 1),
        procs_per_node=112,
        scaling_factor=4.0,
        scaling_dimensions=2,
        max_nodes=128
    )
    seq2 = scaler2.generate_sequence()
    scaler2.print_sequence(seq2)
