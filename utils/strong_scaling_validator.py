"""
Strong Scaling Validator

Validates strong scaling configurations against hardware constraints.
Determines which node counts are mathematically valid for a given problem size.
"""

import logging
from typing import List, Tuple, Dict, Optional
from math import ceil

logger = logging.getLogger(__name__)


class StrongScalingValidator:
    """
    Validates strong scaling configuration for mathematical feasibility.
    
    For strong scaling:
    - Problem size (grid) is FIXED
    - Hardware (cores/node) is FIXED
    - Must find which node counts allow valid MPI decomposition
    """
    
    def __init__(self, 
                 procs_per_node: int,
                 grid_dims: Tuple[int, int, int], 
                 scaling_dims: int = 2,
                 initial_procs: Optional[Tuple[int, int, int]] = None):
        """
        Initialize validator.
        
        Args:
            procs_per_node: Cores available per node (hardware constraint)
            grid_dims: Fixed grid dimensions (nx, ny, nz)
            scaling_dims: Number of dimensions to scale (2 or 3)
            initial_procs: Optional initial MPI decomposition (px, py, pz)
        """
        self.nx, self.ny, self.nz = grid_dims
        self.procs_per_node = procs_per_node
        self.scaling_dims = scaling_dims
        self.initial_procs = initial_procs
        
        # Force pz=1 for 2D or if nz=1
        self.force_pz_one = (scaling_dims == 2 or self.nz == 1)
        
        logger.info(f"Initialized StrongScalingValidator:")
        logger.info(f"  Grid: {self.nx} × {self.ny} × {self.nz} (FIXED)")
        logger.info(f"  Procs/node: {self.procs_per_node} (hardware)")
        logger.info(f"  Scaling dims: {scaling_dims}D")
        if initial_procs:
            logger.info(f"  Initial procs: {initial_procs[0]}×{initial_procs[1]}×{initial_procs[2]}")
    
    def validate_node_count(self, num_nodes: int) -> Tuple[bool, Optional[Tuple], str]:
        """
        Check if a specific node count allows valid decomposition.
        
        Args:
            num_nodes: Number of nodes to test
        
        Returns:
            tuple: (is_valid, decomposition, message)
                - is_valid: True if valid decomposition exists
                - decomposition: (px, py, pz) if valid, None otherwise
                - message: Explanation message
        """
        total_ranks = num_nodes * self.procs_per_node
        
        decomp = self._find_valid_decomposition(total_ranks)
        
        if decomp:
            px, py, pz = decomp
            # Verify perfect divisibility
            if (self.nx % px == 0 and self.ny % py == 0 and 
                (self.nz % pz == 0 or pz == 1)):
                msg = f"✓ {num_nodes} nodes: {total_ranks} ranks = {px}×{py}×{pz}"
                return True, decomp, msg
            else:
                msg = (f"⚠ {num_nodes} nodes: {total_ranks} ranks = {px}×{py}×{pz} "
                      f"(imperfect divisibility)")
                return False, decomp, msg
        else:
            msg = f"❌ {num_nodes} nodes: {total_ranks} ranks - NO VALID DECOMPOSITION"
            return False, None, msg
    
    def validate_sequence(self, max_nodes: int, 
                         sequence_type: str = 'power_of_2') -> Dict:
        """
        Validate an entire node count sequence.
        
        Args:
            max_nodes: Maximum number of nodes
            sequence_type: 'power_of_2' or 'linear' or 'custom'
        
        Returns:
            dict: Validation results with valid/invalid node counts
        """
        if sequence_type == 'power_of_2':
            node_counts = self._generate_power_of_2_sequence(max_nodes)
        else:
            node_counts = list(range(1, max_nodes + 1))
        
        valid_nodes = []
        invalid_nodes = []
        decompositions = {}
        messages = []
        
        logger.info(f"Validating {len(node_counts)} node counts...")
        
        for n in node_counts:
            is_valid, decomp, msg = self.validate_node_count(n)
            messages.append(msg)
            
            if is_valid:
                valid_nodes.append(n)
                decompositions[n] = decomp
                logger.info(f"  {msg}")
            else:
                invalid_nodes.append(n)
                logger.warning(f"  {msg}")
        
        return {
            'valid_nodes': valid_nodes,
            'invalid_nodes': invalid_nodes,
            'decompositions': decompositions,
            'messages': messages,
            'total_tested': len(node_counts),
            'success_rate': len(valid_nodes) / len(node_counts) if node_counts else 0
        }
    
    def validate_node_sequence(self, node_sequence: List[int], 
                               auto_filter: bool = True) -> Tuple[List[int], List[int], Dict]:
        """
        Validate a pre-generated node sequence (called by orchestrator).
        
        Args:
            node_sequence: List of node counts to validate
            auto_filter: If True, filter out invalid nodes automatically
        
        Returns:
            tuple: (valid_nodes, invalid_nodes, details_dict)
        """
        valid_nodes = []
        invalid_nodes = []
        decompositions = {}
        alternatives = {}
        
        logger.info(f"Validating {len(node_sequence)} node counts...")
        
        for n in node_sequence:
            is_valid, decomp, msg = self.validate_node_count(n)
            
            if is_valid:
                valid_nodes.append(n)
                decompositions[n] = decomp
                logger.info(f"  {msg}")
            else:
                invalid_nodes.append(n)
                logger.warning(f"  {msg}")
                
                # Find alternatives for invalid nodes
                alts = self.suggest_alternatives(n, tolerance=0.10)
                if alts:
                    alternatives[n] = alts
        
        details = {
            'decompositions': decompositions,
            'alternatives': alternatives,
            'total_tested': len(node_sequence),
            'success_rate': len(valid_nodes) / len(node_sequence) if node_sequence else 0
        }
        
        return valid_nodes, invalid_nodes, details
    
    def suggest_alternatives(self, target_nodes: int, tolerance: float = 0.2) -> List[Tuple]:
        """
        Suggest alternative node counts near target that would work.
        
        Args:
            target_nodes: Desired number of nodes
            tolerance: Search range (±20% by default)
        
        Returns:
            list: List of (node_count, ranks, decomposition) tuples
        """
        target_ranks = target_nodes * self.procs_per_node
        min_ranks = int(target_ranks * (1 - tolerance))
        max_ranks = int(target_ranks * (1 + tolerance))
        
        alternatives = []
        
        # Find all rank counts in range that have valid decompositions
        for ranks in range(min_ranks, max_ranks + 1):
            decomp = self._find_valid_decomposition(ranks)
            if decomp:
                px, py, pz = decomp
                # Check perfect divisibility
                if (self.nx % px == 0 and self.ny % py == 0 and 
                    (self.nz % pz == 0 or pz == 1)):
                    nodes_needed = ceil(ranks / self.procs_per_node)
                    diff = ranks - target_ranks
                    diff_pct = (diff / target_ranks) * 100
                    alternatives.append((nodes_needed, ranks, decomp, diff, diff_pct))
        
        # Sort by proximity to target
        alternatives.sort(key=lambda x: abs(x[3]))
        
        return alternatives
    
    def suggest_compatible_configs(self, node_sequence: List[int]):
        """
        Print suggestions for incompatible node counts (called by orchestrator).
        
        Args:
            node_sequence: List of node counts that were tested
        """
        logger.info("")
        logger.info("═" * 70)
        logger.info("ALTERNATIVE CONFIGURATIONS FOR INVALID NODE COUNTS")
        logger.info("═" * 70)
        logger.info("")
        
        for n in node_sequence:
            is_valid, _, _ = self.validate_node_count(n)
            
            if not is_valid:
                logger.info(f"Node count {n} (target: {n * self.procs_per_node} ranks):")
                
                alternatives = self.suggest_alternatives(n, tolerance=0.10)
                
                if alternatives:
                    logger.info(f"  Valid alternatives:")
                    for nodes, ranks, decomp, diff, diff_pct in alternatives[:3]:
                        px, py, pz = decomp
                        logger.info(f"    • {nodes} nodes ({ranks} ranks = {px}×{py}×{pz}, {diff_pct:+.1f}%)")
                else:
                    logger.info(f"  No valid alternatives found within ±10%")
                
                logger.info("")
        
        logger.info("═" * 70)
        logger.info("")
    
    def print_validation_report(self, validation_results: Dict):
        """Print comprehensive validation report."""
        print()
        print("="*70)
        print("Strong Scaling Validation Report")
        print("="*70)
        print()
        print(f"Problem Configuration (FIXED):")
        print(f"  Grid dimensions: {self.nx} × {self.ny} × {self.nz}")
        print(f"  Total grid cells: {self.nx * self.ny * self.nz:,}")
        print()
        print(f"Hardware Configuration (FIXED):")
        print(f"  Cores per node: {self.procs_per_node}")
        print()
        print(f"Validation Results:")
        print(f"  Node counts tested: {validation_results['total_tested']}")
        print(f"  Valid configurations: {len(validation_results['valid_nodes'])}")
        print(f"  Invalid configurations: {len(validation_results['invalid_nodes'])}")
        print(f"  Success rate: {validation_results['success_rate']*100:.1f}%")
        print()
        
        if validation_results['valid_nodes']:
            print(f"✓ VALID Node Counts:")
            for n in validation_results['valid_nodes']:
                decomp = validation_results['decompositions'][n]
                px, py, pz = decomp
                ranks = n * self.procs_per_node
                cells_per_proc = (self.nx * self.ny * self.nz) / ranks
                print(f"    {n:3d} nodes: {ranks:5d} ranks = {px:3d}×{py:3d}×{pz:1d} "
                      f"({cells_per_proc:.1f} cells/proc)")
        
        if validation_results['invalid_nodes']:
            print()
            print(f"❌ INVALID Node Counts (will be skipped):")
            for n in validation_results['invalid_nodes']:
                ranks = n * self.procs_per_node
                print(f"    {n:3d} nodes: {ranks:5d} ranks - no valid decomposition")
                
                # Suggest alternatives
                alternatives = self.suggest_alternatives(n, tolerance=0.15)
                if alternatives:
                    print(f"      Alternatives:")
                    for alt_nodes, alt_ranks, alt_decomp, diff, diff_pct in alternatives[:3]:
                        px, py, pz = alt_decomp
                        print(f"        • {alt_nodes} nodes ({alt_ranks} ranks = {px}×{py}×{pz}, "
                              f"{diff_pct:+.1f}%)")
        
        print()
        print("="*70)
        print()
    
    def _find_valid_decomposition(self, total_ranks: int) -> Optional[Tuple[int, int, int]]:
        """
        Find valid MPI decomposition for given rank count.
        
        For 2D (or pz=1):
            total_ranks = px × py × 1
        For 3D:
            total_ranks = px × py × pz
        
        Returns:
            tuple: (px, py, pz) or None if no valid decomposition
        """
        if self.force_pz_one:
            # 2D decomposition
            decomp_2d = self._find_2d_decomposition(total_ranks)
            if decomp_2d:
                return (decomp_2d[0], decomp_2d[1], 1)
            return None
        else:
            # 3D decomposition
            return self._find_3d_decomposition(total_ranks)
    
    def _find_2d_decomposition(self, ranks: int) -> Optional[Tuple[int, int]]:
        """Find px, py such that px×py=ranks, px|nx, py|ny."""
        rank_divisors = self._get_divisors(ranks)
        nx_divisors = set(self._get_divisors(self.nx))
        ny_divisors = set(self._get_divisors(self.ny))
        
        best = None
        best_score = float('inf')
        
        for px in rank_divisors:
            py = ranks // px
            
            if px * py != ranks:
                continue
            
            x_div = (px in nx_divisors)
            y_div = (py in ny_divisors)
            
            if x_div and y_div:
                # Perfect divisibility - calculate aspect ratio score
                grid_aspect = self.nx / self.ny if self.ny > 0 else float('inf')
                decomp_aspect = px / py if py > 0 else float('inf')
                score = abs(grid_aspect - decomp_aspect)
                
                if score < best_score:
                    best_score = score
                    best = (px, py)
        
        return best
    
    def _find_3d_decomposition(self, ranks: int) -> Optional[Tuple[int, int, int]]:
        """Find px, py, pz such that px×py×pz=ranks, divisibility constraints."""
        # TODO: Implement 3D decomposition if needed
        # For now, fall back to 2D with pz=1
        decomp_2d = self._find_2d_decomposition(ranks)
        if decomp_2d:
            return (decomp_2d[0], decomp_2d[1], 1)
        return None
    
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
    
    def _generate_power_of_2_sequence(self, max_nodes: int) -> List[int]:
        """Generate power-of-2 sequence: 1, 2, 4, 8, ..., max_nodes."""
        sequence = []
        n = 1
        while n <= max_nodes:
            sequence.append(n)
            n *= 2
        return sequence
    
    def generate_power_of_2_sequence(self, max_nodes: int) -> List[int]:
        """
        Generate power-of-2 node sequence for strong scaling.
        
        Public method that wraps _generate_power_of_2_sequence.
        
        Args:
            max_nodes: Maximum number of nodes
            
        Returns:
            list: Power-of-2 sequence [1, 2, 4, 8, ..., max_nodes]
        """
        return self._generate_power_of_2_sequence(max_nodes)


if __name__ == "__main__":
    # Test validator
    logging.basicConfig(level=logging.INFO)
    
    print("Testing Strong Scaling Validator")
    print()
    
    # Test case 1: Original problem (840×480, 112 procs/node)
    print("Test 1: Original Configuration (840×480 grid, 112 procs/node)")
    validator1 = StrongScalingValidator(
        procs_per_node=112,
        grid_dims=(840, 480, 1),
        scaling_dims=2
    )
    results1 = validator1.validate_sequence(max_nodes=128)
    validator1.print_validation_report(results1)
    
    # Test case 2: Power-of-2 grid (1024×512, 128 procs/node)
    print("\n" + "="*70 + "\n")
    print("Test 2: Power-of-2 Configuration (1024×512 grid, 128 procs/node)")
    validator2 = StrongScalingValidator(
        procs_per_node=128,
        grid_dims=(1024, 512, 1),
        scaling_dims=2
    )
    results2 = validator2.validate_sequence(max_nodes=128)
    validator2.print_validation_report(results2)
