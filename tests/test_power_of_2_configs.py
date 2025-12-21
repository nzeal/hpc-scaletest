#!/usr/bin/env python3
"""
Verify Power-of-2 Strong Scaling Configurations

Shows which configurations work for complete power-of-2 sequences (1→128 nodes).
"""

import sys

def get_divisors(n):
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

def check_decomposition(ranks, nx, ny):
    """Check if ranks can be decomposed for grid (nx, ny)."""
    nx_divs = set(get_divisors(nx))
    ny_divs = set(get_divisors(ny))
    rank_divs = get_divisors(ranks)
    
    # Find best decomposition
    best = None
    for px in rank_divs:
        py = ranks // px
        if px * py == ranks:
            if px in nx_divs and py in ny_divs:
                return (px, py, True)  # Perfect
            if best is None:
                best = (px, py, False)  # Imperfect
    
    return best if best else (1, 1, False)

def test_config(name, procs_per_node, nx, ny, max_nodes=128):
    """Test a configuration across power-of-2 node counts."""
    print("="*70)
    print(f"Configuration: {name}")
    print("="*70)
    print(f"  procs_per_node: {procs_per_node}")
    print(f"  Grid: {nx} × {ny} × 1")
    print()
    
    nodes = [2**i for i in range(8)]  # 1, 2, 4, 8, 16, 32, 64, 128
    nodes = [n for n in nodes if n <= max_nodes]
    
    all_perfect = True
    
    for n in nodes:
        ranks = n * procs_per_node
        px, py, perfect = check_decomposition(ranks, nx, ny)
        
        status = "✓" if perfect else "❌"
        if not perfect:
            all_perfect = False
        
        cells_x = nx / px if px > 0 else 0
        cells_y = ny / py if py > 0 else 0
        
        print(f"  {status} {n:3d} nodes: {ranks:5d} ranks = {px:4d}×{py:4d} "
              f"(cells/proc: {cells_x:.1f}×{cells_y:.1f})")
    
    print()
    if all_perfect:
        print(f"✅ SUCCESS: All power-of-2 node counts work!")
    else:
        print(f"❌ FAILURE: Some decompositions don't divide evenly")
    print()
    
    return all_perfect

print("\n" + "="*70)
print("Testing Strong Scaling Configurations for Power-of-2 Sequences")
print("="*70)
print()

# Test 1: User's ORIGINAL configuration (BROKEN)
print("Test 1: YOUR ORIGINAL CONFIG (BROKEN)")
test1 = test_config(
    "Original (840×480, 112 procs/node)",
    procs_per_node=112,
    nx=840,
    ny=480
)

# Test 2: RECOMMENDED - Power-of-2 everything
print("Test 2: RECOMMENDED CONFIG (128 procs/node, 1024×512 grid)")
test2 = test_config(
    "Power-of-2 (1024×512, 128 procs/node)",
    procs_per_node=128,
    nx=1024,
    ny=512
)

# Test 3: Alternative - Adjusted for 112 procs/node
print("Test 3: ADJUSTED FOR 112 PROCS/NODE (896×512 grid)")
test3 = test_config(
    "Adjusted (896×512, 112 procs/node)",
    procs_per_node=112,
    nx=896,
    ny=512
)

# Test 4: Alternative - 64 procs/node
print("Test 4: ALTERNATIVE (64 procs/node, 512×512 grid)")
test4 = test_config(
    "Square (512×512, 64 procs/node)",
    procs_per_node=64,
    nx=512,
    ny=512
)

# Summary
print("="*70)
print("SUMMARY")
print("="*70)
print()

configs = [
    ("Original (840×480, 112 procs/node)", test1),
    ("Power-of-2 (1024×512, 128 procs/node)", test2),
    ("Adjusted (896×512, 112 procs/node)", test3),
    ("Square (512×512, 64 procs/node)", test4),
]

for name, result in configs:
    status = "✅ WORKS" if result else "❌ BROKEN"
    print(f"{status}: {name}")

print()
print("RECOMMENDATION:")
print("  Use 'run_strong_power_of_2.yaml' with:")
print("    - procs_per_node: 128")
print("    - initial_cells: [1024, 512, 1]")
print("    - initial_procs: [16, 8, 1]")
print()
print("  This gives PERFECT power-of-2 scaling from 1 to 128 nodes!")
print()
