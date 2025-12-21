#!/usr/bin/env python3
"""
HPC Strong Scaling Configuration Validator

Smart tool that:
1. Auto-detects hardware (cores/node)
2. Validates your FIXED problem size
3. Shows which node counts work
4. Suggests alternatives if needed

Usage:
    python validate_strong_scaling.py --grid 840 480 1 --nodes 128
    python validate_strong_scaling.py --config run_strong.yaml
"""

import sys
import argparse
import logging

# Add parent to path
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.hardware_detector import HardwareDetector
from utils.strong_scaling_validator import StrongScalingValidator

logging.basicConfig(level=logging.WARNING)  # Only show warnings/errors


def parse_args():
    parser = argparse.ArgumentParser(
        description='Validate strong scaling configuration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate with auto-detected hardware
  python validate_strong_scaling.py --grid 840 480 1 --nodes 128
  
  # Specify hardware explicitly
  python validate_strong_scaling.py --grid 840 480 1 --nodes 128 --procs-per-node 112
  
  # Read from YAML config
  python validate_strong_scaling.py --config run_strong.yaml
        """
    )
    
    parser.add_argument('--grid', nargs=3, type=int, metavar=('NX', 'NY', 'NZ'),
                       help='Grid dimensions (nx ny nz)')
    parser.add_argument('--nodes', type=int, metavar='N',
                       help='Maximum number of nodes')
    parser.add_argument('--procs-per-node', type=int, metavar='P',
                       help='Processes per node (auto-detect if not specified)')
    parser.add_argument('--scaling-dims', type=int, choices=[2, 3], default=2,
                       help='Scaling dimensions (default: 2)')
    parser.add_argument('--config', type=str, metavar='FILE',
                       help='Read configuration from YAML file')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    return parser.parse_args()


def read_yaml_config(config_file):
    """Read configuration from YAML file."""
    try:
        import yaml
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        grid = config.get('initial_cells', [0, 0, 0])
        nodes = config.get('nodes', 1)
        procs_per_node = config.get('procs_per_node', None)
        scaling_dims = config.get('scaling_dimensions', 2)
        
        return grid, nodes, procs_per_node, scaling_dims
    except ImportError:
        print("❌ Error: PyYAML not installed. Install with: pip install pyyaml")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error reading config file: {e}")
        sys.exit(1)


def main():
    args = parse_args()
    
    if args.verbose:
        logging.basicConfig(level=logging.INFO, force=True)
    
    # Get configuration
    if args.config:
        grid, max_nodes, procs_per_node, scaling_dims = read_yaml_config(args.config)
        print(f"Configuration read from: {args.config}")
    else:
        if not args.grid or not args.nodes:
            print("❌ Error: Must specify either --config or both --grid and --nodes")
            print("Run with --help for usage information")
            sys.exit(1)
        
        grid = tuple(args.grid)
        max_nodes = args.nodes
        procs_per_node = args.procs_per_node
        scaling_dims = args.scaling_dims
    
    print()
    print("="*70)
    print("HPC Strong Scaling Configuration Validator")
    print("="*70)
    print()
    
    # Step 1: Detect/confirm hardware
    if procs_per_node is None:
        print("Step 1: Auto-detecting hardware...")
        print("-"*70)
        detector = HardwareDetector()
        hw_config = detector.detect_all()
        
        if hw_config['cores_per_node']:
            procs_per_node = hw_config['cores_per_node']
            print(f"✓ Detected: {procs_per_node} cores/node")
        else:
            print("❌ Could not auto-detect cores per node")
            print("Please specify manually with --procs-per-node")
            sys.exit(1)
    else:
        print("Step 1: Hardware Configuration")
        print("-"*70)
        print(f"✓ Using specified: {procs_per_node} cores/node")
    
    print()
    
    # Step 2: Problem configuration
    print("Step 2: Problem Configuration (FIXED for strong scaling)")
    print("-"*70)
    nx, ny, nz = grid
    print(f"  Grid dimensions: {nx} × {ny} × {nz}")
    print(f"  Total cells: {nx * ny * nz:,}")
    print(f"  Scaling mode: {scaling_dims}D")
    print()
    
    # Step 3: Validate
    print("Step 3: Validating Node Count Sequence")
    print("-"*70)
    print(f"  Testing: 1 → 2 → 4 → 8 → ... → {max_nodes} nodes")
    print()
    
    validator = StrongScalingValidator(
        procs_per_node=procs_per_node,
        grid_dims=grid,
        scaling_dims=scaling_dims
    )
    
    results = validator.validate_sequence(max_nodes=max_nodes)
    
    # Step 4: Report
    validator.print_validation_report(results)
    
    # Step 5: Recommendations
    if results['invalid_nodes']:
        print("="*70)
        print("RECOMMENDATIONS")
        print("="*70)
        print()
        print("Your configuration has invalid node counts. You have 3 options:")
        print()
        print("OPTION 1: Use only VALID node counts (RECOMMENDED)")
        print(f"  Update your YAML with:")
        print(f"    node_counts: {results['valid_nodes']}")
        print()
        print("OPTION 2: Adjust problem size (changes your science!)")
        print(f"  Find grid dimensions compatible with your hardware")
        print(f"  Run: python find_compatible_grid.py --procs-per-node {procs_per_node}")
        print()
        print("OPTION 3: Framework auto-skip (less control)")
        print(f"  Framework will automatically skip invalid node counts")
        print(f"  and log warnings")
        print()
        
        # Calculate success rate
        if results['success_rate'] < 0.5:
            print("⚠ WARNING: Less than 50% of node counts are valid!")
            print("  Strongly recommend OPTION 1 or OPTION 2")
        
        sys.exit(1)
    else:
        print("="*70)
        print("✅ SUCCESS!")
        print("="*70)
        print()
        print("All node counts in your sequence have valid decompositions.")
        print("Your configuration is ready to run!")
        print()
        sys.exit(0)


if __name__ == "__main__":
    main()
