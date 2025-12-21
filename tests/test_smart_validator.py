#!/usr/bin/env python3
"""
Test Smart Strong Scaling Validator

Demonstrates the intelligent pre-flight validation system that:
1. Detects hardware constraints (procs_per_node)
2. Analyzes problem size (grid dimensions)
3. Predicts which node counts will work/fail
4. Auto-filters invalid configurations
5. Suggests compatible alternatives
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.strong_scaling_validator import StrongScalingValidator

print("="*70)
print("INTELLIGENT STRONG SCALING VALIDATOR TEST")
print("="*70)
print()

# Test Case 1: User's ORIGINAL configuration (some nodes fail)
print("Test Case 1: YOUR ORIGINAL CONFIGURATION")
print("-"*70)
print()

validator1 = StrongScalingValidator(
    procs_per_node=112,  # Leonardo Booster (FIXED by hardware)
    grid_dims=(840, 480, 1),  # Your grid (FIXED for strong scaling)
    scaling_dims=2
)

node_sequence = [1, 2, 4, 8, 16, 32, 64, 128]
valid, invalid, details = validator1.validate_node_sequence(node_sequence, auto_filter=True)

print()
print(f"Result: {len(valid)}/{len(node_sequence)} node counts are valid")
print(f"  Valid:   {valid}")
print(f"  Invalid: {invalid}")
print()

if invalid:
    print("⚠ PROBLEM DETECTED: Some node counts will fail!")
    print("  Framework will AUTO-FILTER and skip these nodes")
    print()
    validator1.suggest_compatible_configs(node_sequence)

input("Press ENTER to continue to next test case...")
print()
print()

# Test Case 2: Power-of-2 configuration (all nodes work)
print("Test Case 2: POWER-OF-2 CONFIGURATION")
print("-"*70)
print()

validator2 = StrongScalingValidator(
    procs_per_node=128,  # Power of 2
    grid_dims=(1024, 512, 1),  # Powers of 2
    scaling_dims=2
)

valid2, invalid2, details2 = validator2.validate_node_sequence(node_sequence, auto_filter=False)

print()
print(f"Result: {len(valid2)}/{len(node_sequence)} node counts are valid")
if invalid2:
    print(f"  Invalid: {invalid2}")
else:
    print("  ✅ ALL NODE COUNTS WORK!")
print()

input("Press ENTER to continue to next test case...")
print()
print()

# Test Case 3: Show what happens WITHOUT validator (old behavior)
print("Test Case 3: WITHOUT VALIDATOR (OLD BEHAVIOR)")
print("-"*70)
print()
print("OLD FRAMEWORK BEHAVIOR:")
print("  1. Generate jobs for ALL node counts: 1, 2, 4, 8, 16, 32, 64, 128")
print("  2. Submit ALL jobs (no filtering)")
print("  3. Jobs for 32, 64, 128 nodes CRASH at runtime with:")
print("     'MPI_Cart_create error: communicator (3584) < topology (3600)'")
print("  4. User discovers failure AFTER wasting time/resources")
print()
print("NEW FRAMEWORK BEHAVIOR (v1.4.0+):")
print("  1. PRE-FLIGHT validation BEFORE generating any jobs")
print("  2. Detect: 32, 64, 128 will fail")
print("  3. AUTO-FILTER: Only generate jobs for 1, 2, 4, 8, 16")
print("  4. Suggest alternatives (different grid dimensions)")
print("  5. User sees problems BEFORE wasting resources!")
print()

input("Press ENTER to see demonstration...")
print()
print()

# Test Case 4: Detection for a specific problematic node count
print("Test Case 4: DETAILED ANALYSIS OF 32 NODES")
print("-"*70)
print()

validator4 = StrongScalingValidator(
    procs_per_node=112,
    grid_dims=(840, 480, 1),
    scaling_dims=2
)

total_ranks = 32 * 112
is_valid, decomp, reason = validator4.can_decompose(total_ranks)

print(f"Analyzing: 32 nodes × 112 procs/node = {total_ranks} ranks")
print(f"Grid: 840 × 480 × 1 (FIXED for strong scaling)")
print()
print(f"Question: Can {total_ranks} ranks be decomposed as px × py × 1")
print(f"          where px divides 840 and py divides 480?")
print()

if is_valid:
    print(f"✓ YES: {decomp}")
else:
    print(f"❌ NO: {reason}")
    print()
    print("Mathematical explanation:")
    print(f"  {total_ranks} = 2^9 × 7 = 512 × 7")
    print(f"  840 divisors: {validator4._get_divisors(840)[:10]}...")
    print(f"  480 divisors: {validator4._get_divisors(480)[:10]}...")
    print()
    print("  No pair (px, py) exists where:")
    print(f"    - px × py = {total_ranks}")
    print("    - px divides 840")
    print("    - py divides 480")
    print()
    print("  Closest match: 60 × 60 = 3600")
    print("    Both 60 divides 840 and 60 divides 480")
    print("    But 3600 ≠ 3584 → MPI_Cart_create error!")

print()
input("Press ENTER to see summary...")
print()
print()

# Summary
print("="*70)
print("SUMMARY: INTELLIGENT VALIDATION")
print("="*70)
print()
print("KEY FEATURES:")
print()
print("1. HARDWARE-AWARE")
print("   • Knows procs_per_node is FIXED (can't change)")
print("   • Leonardo Booster: 112 cores/node")
print()
print("2. PROBLEM-AWARE")
print("   • Knows grid dimensions are FIXED (strong scaling)")
print("   • Your grid: 840 × 480 × 1")
print()
print("3. PREDICTIVE")
print("   • Calculates ALL valid/invalid node counts UPFRONT")
print("   • No runtime surprises!")
print()
print("4. AUTO-FILTERING")
print("   • Automatically skips invalid node counts")
print("   • Only generates jobs that will work")
print()
print("5. HELPFUL SUGGESTIONS")
print("   • Suggests alternative configurations")
print("   • Shows compatible grid dimensions")
print()
print("BENEFITS:")
print("  ✓ No wasted compute time on failing jobs")
print("  ✓ No MPI_Cart_create errors at runtime")
print("  ✓ Clear warnings before job submission")
print("  ✓ Actionable suggestions for fixes")
print()
print("="*70)
print()
print("The validator is now integrated into the framework!")
print("It runs automatically for all strong scaling tests.")
print()
