#!/usr/bin/env python3
"""
Comprehensive MPI + GPU Detection Test

Tests both MPI implementation detection and GPU configuration detection,
then shows adapted launch commands.

Usage:
  python3 test_mpi_gpu_detection.py [partition]
  
Example:
  python3 test_mpi_gpu_detection.py boost_usr_prod
"""

import sys
import os
import logging

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import detection modules
try:
    from utils.mpi_detector import MPIDetector
    HAS_MPI_DETECTOR = True
except ImportError:
    HAS_MPI_DETECTOR = False
    print("⚠ Warning: MPI detector not available")

try:
    from utils.advanced_gpu_manager import AdvancedGPUManager
    HAS_GPU_MANAGER = True
except ImportError:
    HAS_GPU_MANAGER = False
    print("⚠ Warning: GPU manager not available")

logging.basicConfig(level=logging.INFO, format='%(message)s')


def print_section(title):
    """Print a section header."""
    print()
    print("=" * 70)
    print(title)
    print("=" * 70)
    print()


def test_mpi_detection():
    """Test MPI implementation detection."""
    if not HAS_MPI_DETECTOR:
        print("❌ MPI detector not available")
        return None
    
    print_section("Step 1: MPI Implementation Detection")
    
    detector = MPIDetector()
    info = detector.detect()
    
    print("MPI Implementation: {}".format(info.implementation.upper()))
    print("MPI Version:        {}".format(info.version))
    if info.variant:
        print("MPI Variant:        {}".format(info.variant))
    print()
    
    print("Supported Features:")
    if info.features:
        for feature in sorted(info.features):
            print("  ✓ {}".format(feature))
    else:
        print("  (none - uses standard MPI only)")
    print()
    
    print("Capabilities:")
    print("  --report-bindings:  {}".format("YES ✓" if info.supports_report_bindings() else "NO"))
    print("  ppr mapping:        {}".format("YES ✓" if info.supports_ppr_mapping() else "NO"))
    print("  Preferred launcher: {}".format(info.get_launcher_name()))
    
    return info


def test_gpu_detection(partition):
    """Test GPU configuration detection."""
    if not HAS_GPU_MANAGER:
        print("❌ GPU manager not available")
        return None
    
    print_section("Step 2: GPU Hardware Detection")
    
    manager = AdvancedGPUManager()
    config = manager.detect_gpu_node_config(partition)
    
    if not config:
        print("❌ Could not detect GPU configuration for partition '{}'".format(partition))
        return None
    
    print("Partition:          {}".format(config.partition))
    print("GPUs per node:      {}".format(config.gpus_per_node))
    print("GPU model:          {}".format(config.gpu_model))
    print("CPUs per node:      {}".format(config.cpus_per_node))
    print("Cores per GPU:      {}".format(config.cores_per_gpu))
    print("Memory:             {} GB/node".format(config.memory_per_node_gb))
    
    return config


def show_adapted_commands(mpi_info, gpu_config):
    """Show MPI-adapted launch commands for detected GPU configuration."""
    if not gpu_config or not HAS_GPU_MANAGER:
        return
    
    print_section("Step 3: MPI-Adapted Launch Commands")
    
    manager = AdvancedGPUManager()
    manager.node_config = gpu_config
    
    print("Target Configuration:")
    print("  {} GPUs per node".format(gpu_config.gpus_per_node))
    print("  {} cores per GPU".format(gpu_config.cores_per_gpu))
    print("  {} MPI implementation".format(
        mpi_info.implementation.upper() if mpi_info else "unknown"
    ))
    print()
    
    # Production command (no verbose)
    print("Production Command (recommended):")
    cmd_prod = manager.generate_mpirun_command(
        executable="./app",
        args="input.txt",
        verbose=False
    )
    print("  {}".format(cmd_prod))
    print()
    
    # Debug command (verbose if supported)
    print("Debug Command (with verbose output, if supported):")
    cmd_debug = manager.generate_mpirun_command(
        executable="./app",
        args="input.txt",
        verbose=True
    )
    print("  {}".format(cmd_debug))
    print()
    
    # Explain differences
    if mpi_info:
        print("MPI-Specific Adaptations:")
        if mpi_info.implementation == 'openmpi':
            print("  ✓ Using ppr mapping: --map-by ppr:{}:node --bind-to core --cpus-per-proc {}".format(
                gpu_config.gpus_per_node, gpu_config.cores_per_gpu
            ))
            print("  ✓ Verbose mode adds: --report-bindings")
        elif mpi_info.implementation == 'intelmpi':
            print("  ✓ Using standard MPI syntax (no ppr mapping)")
            print("  ✓ Process placement handled by SLURM")
            print("  ℹ --report-bindings not supported (skipped)")
        else:
            print("  ✓ Using MPI-agnostic syntax")
            print("  ✓ Process placement handled by SLURM")


def show_slurm_directives(gpu_config):
    """Show SLURM directives."""
    if not gpu_config or not HAS_GPU_MANAGER:
        return
    
    print_section("Step 4: SLURM Job Directives")
    
    manager = AdvancedGPUManager()
    manager.node_config = gpu_config
    
    directives = manager.generate_slurm_directives(num_nodes=1)
    
    print("For 1 node:")
    print("  #SBATCH --nodes=1")
    print("  #SBATCH --partition={}".format(directives['partition']))
    print("  #SBATCH --ntasks-per-node={}  # CPU allocation".format(directives['ntasks-per-node']))
    print("  #SBATCH --cpus-per-task={}".format(directives['cpus-per-task']))
    print("  #SBATCH --gres={}        # GPU allocation".format(directives['gres']))
    print()
    
    print("Key Point:")
    print("  --ntasks-per-node={} allocates {} CPU cores".format(
        directives['ntasks-per-node'], gpu_config.cpus_per_node
    ))
    print("  Actual MPI tasks: {} (1 per GPU)".format(gpu_config.gpus_per_node))
    print("  Each task gets: {} cores".format(gpu_config.cores_per_gpu))


def show_summary(mpi_info, gpu_config):
    """Show summary and recommendations."""
    print_section("Summary & Recommendations")
    
    print("Detection Status:")
    print("  MPI:  {}".format("✓ Detected" if mpi_info else "❌ Not detected"))
    print("  GPU:  {}".format("✓ Detected" if gpu_config else "❌ Not detected"))
    print()
    
    if mpi_info and gpu_config:
        print("Your System Configuration:")
        print("  MPI:  {} {}".format(mpi_info.implementation.upper(), mpi_info.version))
        print("  GPU:  {} × {} ({})".format(
            gpu_config.gpus_per_node,
            gpu_config.gpu_model,
            gpu_config.partition
        ))
        print()
        
        print("Framework Behavior:")
        if mpi_info.implementation == 'openmpi':
            print("  ✓ Will use OpenMPI-specific optimizations")
            print("  ✓ Will add --map-by ppr:N:node --bind-to core for optimal GPU mapping")
            print("  ✓ Can use --report-bindings for debugging")
        elif mpi_info.implementation == 'intelmpi':
            print("  ✓ Will use Intel MPI compatible syntax")
            print("  ✓ Will rely on SLURM for process placement (automatic)")
            print("  ✓ Will skip OpenMPI-specific flags (--report-bindings, ppr mapping)")
        else:
            print("  ✓ Will use MPI-agnostic syntax")
            print("  ✓ Will work with any MPI implementation")
        print()
        
        print("Next Steps:")
        print("  1. Update your config.yaml:")
        print("     partition: {}".format(gpu_config.partition))
        print("     launcher: mpirun  # Framework auto-adapts")
        print()
        print("  2. Run framework:")
        print("     python3 hpc_auto.py --config your_config.yaml --dry-run")
        print()
        print("  3. Check generated job script:")
        print("     cat output/*/node1/job.sh")
        print()
        print("  4. Verify launcher command matches your MPI:")
        if mpi_info.implementation == 'openmpi':
            print("     Should contain: --map-by ppr:{}:node".format(
                gpu_config.gpus_per_node
            ))
        else:
            print("     Should NOT contain: --map-by or --report-bindings")
        print()
        
        print("✓ Your setup is ready!")
    else:
        print("⚠ Detection incomplete. Please check:")
        if not mpi_info:
            print("  • Is MPI module loaded? (module load intelmpi/openmpi)")
        if not gpu_config:
            print("  • Is partition name correct? (sinfo -o '%R')")


def main():
    """Main test function."""
    if len(sys.argv) < 2:
        print("Usage: python3 test_mpi_gpu_detection.py <partition>")
        print()
        print("Examples:")
        print("  python3 test_mpi_gpu_detection.py boost_usr_prod")
        print("  python3 test_mpi_gpu_detection.py booster")
        print("  python3 test_mpi_gpu_detection.py dcgp")
        print()
        sys.exit(1)
    
    partition = sys.argv[1]
    
    print()
    print("=" * 70)
    print("HPC-ScaleTest: MPI + GPU Detection Test")
    print("=" * 70)
    print()
    print("Testing partition: '{}'".format(partition))
    
    # Test MPI detection
    mpi_info = test_mpi_detection()
    
    # Test GPU detection
    gpu_config = test_gpu_detection(partition)
    
    # Show adapted commands
    if gpu_config:
        show_adapted_commands(mpi_info, gpu_config)
        show_slurm_directives(gpu_config)
    
    # Show summary
    show_summary(mpi_info, gpu_config)
    
    print()
    print("=" * 70)
    print()


if __name__ == '__main__':
    main()
