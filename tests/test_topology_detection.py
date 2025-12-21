#!/usr/bin/env python3
"""
Topology Detection and MPI Command Generation Test

This script demonstrates and tests the centralized topology detection
and MPI command generation capabilities of HPC-ScaleTest v2.0.0.

It verifies:
1. Topology detection from multiple sources
2. MPI mapping computation
3. Correct mpirun command syntax
4. GPU binding script generation

Usage:
    python3 tests/test_topology_detection.py [--partition PARTITION]

Example output for a 32-core, 4-GPU node:
    Topology: 32 cores, 4 GPUs (nvidia)
    Mapping: 4 ranks/node, 8 cores/rank, 1 GPU/rank
    mpirun -np 4 --map-by ppr:4:node --bind-to core --cpus-per-proc 8 ./app
"""

import os
import sys
import argparse
import logging

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the centralized modules
try:
    from core.topology import (
        TopologyDetector, NodeTopology, MPIMapping,
        get_topology_detector, GPUVendor, DetectionContext
    )
    HAS_TOPOLOGY = True
except ImportError as e:
    print(f"WARNING: Could not import core.topology: {e}")
    HAS_TOPOLOGY = False

try:
    from core.mpi_command import (
        MPICommandGenerator, MPIDetector, MPIInfo,
        generate_gpu_binding_script, MPIImplementation
    )
    HAS_MPI_COMMAND = True
except ImportError as e:
    print(f"WARNING: Could not import core.mpi_command: {e}")
    HAS_MPI_COMMAND = False


def test_topology_detection(partition=None):
    """Test topology detection."""
    print("\n" + "=" * 70)
    print(" TOPOLOGY DETECTION TEST")
    print("=" * 70)
    
    if not HAS_TOPOLOGY:
        print("SKIP: core.topology module not available")
        return None
    
    detector = TopologyDetector()
    
    # Show detection context
    print(f"\nDetection context: {detector.context.value}")
    
    # Check for Slurm environment
    slurm_vars = [
        'SLURM_JOB_ID', 'SLURM_CPUS_ON_NODE', 'SLURM_GPUS_ON_NODE',
        'SLURM_JOB_PARTITION', 'SLURM_NTASKS_PER_NODE'
    ]
    print("\nSlurm environment variables:")
    for var in slurm_vars:
        value = os.environ.get(var, "(not set)")
        print(f"  {var}: {value}")
    
    # Detect topology
    print(f"\nDetecting topology for partition: {partition or '(auto)'}")
    try:
        topology = detector.detect(partition)
        print(f"\n✓ Detection successful!")
        print(f"  Detection method: {topology.detection_method}")
        print(f"  CPU cores: {topology.cpu_cores}")
        print(f"  Physical cores: {topology.physical_cores}")
        print(f"  Threads per core: {topology.threads_per_core}")
        print(f"  GPUs: {topology.gpus}")
        print(f"  GPU vendor: {topology.gpu_vendor.value}")
        print(f"  GPU model: {topology.gpu_model or '(unknown)'}")
        return topology
    except Exception as e:
        print(f"\n✗ Detection failed: {e}")
        return None


def test_mpi_mapping(topology, num_nodes=2):
    """Test MPI mapping computation."""
    print("\n" + "=" * 70)
    print(" MPI MAPPING TEST")
    print("=" * 70)
    
    if not HAS_TOPOLOGY:
        print("SKIP: core.topology module not available")
        return None
    
    if topology is None:
        print("SKIP: No topology available")
        return None
    
    detector = get_topology_detector()
    
    print(f"\nComputing MPI mapping for {num_nodes} nodes...")
    try:
        mapping = detector.compute_mpi_mapping(topology, num_nodes=num_nodes)
        print(f"\n✓ Mapping computed!")
        print(f"  Total ranks: {mapping.total_ranks}")
        print(f"  Ranks per node: {mapping.ranks_per_node}")
        print(f"  Cores per rank: {mapping.cores_per_rank}")
        print(f"  GPUs per rank: {mapping.gpus_per_rank}")
        
        # Validation
        is_valid, error = mapping.validate(topology)
        if is_valid:
            print(f"  Validation: ✓ PASSED")
        else:
            print(f"  Validation: ✗ FAILED - {error}")
        
        return mapping
    except Exception as e:
        print(f"\n✗ Mapping computation failed: {e}")
        return None


def test_mpi_command_generation(topology, mapping):
    """Test MPI command generation."""
    print("\n" + "=" * 70)
    print(" MPI COMMAND GENERATION TEST")
    print("=" * 70)
    
    if not HAS_MPI_COMMAND:
        print("SKIP: core.mpi_command module not available")
        return
    
    if topology is None or mapping is None:
        print("SKIP: No topology/mapping available")
        return
    
    # Detect MPI
    print("\nDetecting MPI implementation...")
    mpi_detector = MPIDetector()
    mpi_info = mpi_detector.detect()
    print(f"  Implementation: {mpi_info.implementation.value}")
    print(f"  Version: {mpi_info.version or '(unknown)'}")
    print(f"  Launcher: {mpi_info.launcher}")
    print(f"  Supports ppr: {mpi_info.supports_ppr}")
    print(f"  Supports bind-to: {mpi_info.supports_bind_to}")
    print(f"  Supports cpus-per-proc: {mpi_info.supports_cpus_per_proc}")
    
    # Generate command
    print("\nGenerating MPI command...")
    generator = MPICommandGenerator(mpi_info)
    num_nodes = mapping.total_ranks // mapping.ranks_per_node
    cmd = generator.generate(
        topology=topology,
        mapping=mapping,
        executable="./my_application",
        args=["--input", "config.inp"],
        num_nodes=num_nodes,
        verbose=True,
        gpu_binding_script="./bind_gpu.sh" if topology.gpus > 0 else None,
    )
    
    cmd_str = ' '.join(cmd)
    print(f"\nGenerated command:")
    print(f"  {cmd_str}")
    
    # Verify correct syntax
    print("\nSyntax verification:")
    
    # Check for INCORRECT :PE= syntax
    if ':PE=' in cmd_str:
        print("  ✗ ERROR: Found incorrect :PE= syntax!")
        print("    The correct syntax is:")
        print("      --map-by ppr:N:node --bind-to core --cpus-per-proc M")
        print("    NOT:")
        print("      --map-by ppr:N:node:PE=M")
    else:
        print("  ✓ No incorrect :PE= syntax found")
    
    # Check for correct ppr mapping
    if '--map-by' in cmd_str and 'ppr:' in cmd_str:
        print("  ✓ Correct --map-by ppr:N:node syntax")
    
    # Check for correct binding
    if '--bind-to core' in cmd_str:
        print("  ✓ Correct --bind-to core syntax")
    
    # Check for cpus-per-proc
    if '--cpus-per-proc' in cmd_str:
        print("  ✓ Correct --cpus-per-proc syntax")
    
    return cmd


def test_example_configurations():
    """Test with example configurations."""
    print("\n" + "=" * 70)
    print(" EXAMPLE CONFIGURATION TEST (Leonardo Booster-style)")
    print("=" * 70)
    
    if not HAS_TOPOLOGY or not HAS_MPI_COMMAND:
        print("SKIP: Required modules not available")
        return
    
    # Example: Leonardo Booster-like node (32 cores, 4 GPUs)
    print("\nExample: 32 CPU cores, 4 A100 GPUs")
    
    topology = NodeTopology(
        cpu_cores=32,
        gpus=4,
        gpu_vendor=GPUVendor.NVIDIA,
        gpu_model="A100-SXM4-64GB",
        detection_method="example_configuration"
    )
    
    detector = get_topology_detector()
    
    # Create explicit OpenMPI info for testing
    openmpi_info = MPIInfo(
        implementation=MPIImplementation.OPENMPI,
        version="4.1.0",
        launcher="mpirun",
        supports_ppr=True,
        supports_bind_to=True,
        supports_cpus_per_proc=True,
        supports_report_bindings=True,
    )
    generator = MPICommandGenerator(openmpi_info)
    
    # Single node
    print("\nCase 1: Single node job")
    mapping = detector.compute_mpi_mapping(topology, num_nodes=1)
    print(f"  Expected: 4 ranks/node, 8 cores/rank")
    print(f"  Computed: {mapping.ranks_per_node} ranks/node, {mapping.cores_per_rank} cores/rank")
    
    cmd = generator.generate(
        topology=topology,
        mapping=mapping,
        executable="./app",
        num_nodes=1,
    )
    print(f"  Command: {' '.join(cmd)}")
    
    # Verify
    if mapping.ranks_per_node == 4 and mapping.cores_per_rank == 8:
        print("  ✓ PASS: Correct mapping computed")
    else:
        print("  ✗ FAIL: Incorrect mapping")
    
    # Multi-node
    print("\nCase 2: 4-node job")
    mapping = detector.compute_mpi_mapping(topology, num_nodes=4)
    print(f"  Expected: 16 total ranks (4 × 4)")
    print(f"  Computed: {mapping.total_ranks} total ranks")
    
    cmd = generator.generate(
        topology=topology,
        mapping=mapping,
        executable="./app",
        num_nodes=4,
    )
    print(f"  Command: {' '.join(cmd)}")
    
    # Verify
    if mapping.total_ranks == 16:
        print("  ✓ PASS: Correct total ranks")
    else:
        print("  ✗ FAIL: Incorrect total ranks")
    
    # Verify correct mpirun syntax
    cmd_str = ' '.join(cmd)
    print("\nVerifying mpirun syntax (OpenMPI):")
    print(f"  Command: {cmd_str}")
    
    expected_patterns = [
        ('-np 16', 'Correct -np'),
        ('--map-by ppr:4:node', 'Correct ppr mapping'),
        ('--bind-to core', 'Correct binding'),
        ('--cpus-per-proc 8', 'Correct cores per process'),
    ]
    
    all_pass = True
    for pattern, description in expected_patterns:
        if pattern in cmd_str:
            print(f"  ✓ {description}: {pattern}")
        else:
            print(f"  ✗ Missing {description}: {pattern}")
            all_pass = False
    
    if ':PE=' in cmd_str:
        print("  ✗ ERROR: Incorrect :PE= syntax found!")
        all_pass = False
    else:
        print("  ✓ No incorrect :PE= syntax")
    
    return all_pass


def main():
    parser = argparse.ArgumentParser(
        description='Test topology detection and MPI command generation'
    )
    parser.add_argument(
        '--partition', '-p',
        help='Slurm partition to test (optional)',
        default=None
    )
    parser.add_argument(
        '--nodes', '-n',
        type=int,
        default=2,
        help='Number of nodes to test with (default: 2)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    args = parser.parse_args()
    
    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(
        level=log_level,
        format='%(levelname)s - %(name)s - %(message)s'
    )
    
    print("=" * 70)
    print(" HPC-ScaleTest v2.0.0 Topology Detection Test")
    print("=" * 70)
    print(f"\nTest configuration:")
    print(f"  Partition: {args.partition or '(auto-detect)'}")
    print(f"  Nodes: {args.nodes}")
    
    # Run tests
    topology = test_topology_detection(args.partition)
    mapping = test_mpi_mapping(topology, args.nodes)
    test_mpi_command_generation(topology, mapping)
    all_pass = test_example_configurations()
    
    print("\n" + "=" * 70)
    print(" TEST SUMMARY")
    print("=" * 70)
    if all_pass:
        print("\n✓ All example configuration tests PASSED")
    else:
        print("\n✗ Some tests FAILED")
    print()


if __name__ == '__main__':
    main()
