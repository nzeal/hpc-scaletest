#!/usr/bin/env python3
"""
HPC Auto - Automated HPC Scaling Framework CLI
Simple command-line interface for end-to-end HPC scaling tests.
"""

import argparse
import sys
import logging
from pathlib import Path
from typing import Optional

from engine.orchestrator import HPCOrchestrator, OrchestratorConfig
from utils.logging_config import setup_logging
from utils.system_loader import SystemConfigLoader, auto_load_system_config


def create_parser():
    """Create argument parser."""
    parser = argparse.ArgumentParser(
        description="Automated HPC Scaling Framework - Simple, intelligent scaling tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Strong scaling test with auto-detected system config
  python hpc_auto.py /path/to/code --scaling strong --nodes 8

  # Weak scaling test from Git repository
  python hpc_auto.py https://github.com/user/hpc-app.git --scaling weak --nodes 16

  # Use specific system configuration file
  python hpc_auto.py /path/to/code --system-config leonardo_system.py \\
      --partition-name booster --nodes 16

  # Use custom launcher from system config
  python hpc_auto.py /path/to/code --system-config leonardo_system.py \\
      --partition-name booster --use-custom-launcher mpirun-mapby --nodes 8

  # GPU-enabled test with system config
  python hpc_auto.py /path/to/cuda-app --hardware gpu --partition-name booster --nodes 4

  # Manual configuration (disable system config)
  python hpc_auto.py /path/to/code --no-system-config --scaling strong --nodes 32 \\
      --procs-per-node 128 --partition gpu_prod --account project123
        """
    )
    
    # Required arguments
    parser.add_argument(
        'source',
        help='Path to source code or Git repository URL'
    )
    
    # Scaling configuration
    scaling_group = parser.add_argument_group('Scaling Configuration')
    scaling_group.add_argument(
        '--scaling', '-s',
        choices=['strong', 'weak'],
        default='strong',
        help='Scaling type (default: strong)'
    )
    scaling_group.add_argument(
        '--nodes', '-n',
        type=int,
        default=4,
        help='Maximum number of nodes to test (default: 4)'
    )
    scaling_group.add_argument(
        '--initial-procs',
        type=str,
        default='2,2,2',
        help='Initial process decomposition as x,y,z (default: 2,2,2)'
    )
    scaling_group.add_argument(
        '--initial-domain',
        type=str,
        help='Initial domain size as x,y,z (for weak scaling)'
    )
    scaling_group.add_argument(
        '--initial-cells',
        type=str,
        help='Initial cell count as x,y,z (for weak scaling)'
    )
    
    # Hardware configuration
    hardware_group = parser.add_argument_group('Hardware Configuration')
    hardware_group.add_argument(
        '--hardware',
        choices=['cpu', 'gpu'],
        default='cpu',
        help='Target hardware type (default: cpu)'
    )
    hardware_group.add_argument(
        '--procs-per-node',
        type=int,
        default=128,
        help='Processes per node (default: 128)'
    )
    hardware_group.add_argument(
        '--gpus-per-node',
        type=int,
        default=0,
        help='GPUs per node (default: 0, auto-set to 1 for --hardware gpu)'
    )
    
    # Resource configuration
    resource_group = parser.add_argument_group('Resource Configuration')
    resource_group.add_argument(
        '--time-limit',
        default='02:00:00',
        help='Job time limit (default: 02:00:00)'
    )
    resource_group.add_argument(
        '--partition',
        default='X_usr_prod',
        help='Partition/queue name (default: X_usr_prod)'
    )
    resource_group.add_argument(
        '--account',
        default='cin_X',
        help='Account/project name (default: cin_X)'
    )
    
    # System configuration
    system_group = parser.add_argument_group('System Configuration (NEW)')
    system_group.add_argument(
        '--system-config',
        type=Path,
        help='Path to system configuration file (e.g., leonardo_system.py)'
    )
    system_group.add_argument(
        '--partition-name',
        type=str,
        help='Partition name from system config (e.g., booster, dcgp)'
    )
    system_group.add_argument(
        '--use-custom-launcher',
        type=str,
        help='Use custom launcher from system config (e.g., mpirun-mapby, mpirun-nsys)'
    )
    system_group.add_argument(
        '--no-system-config',
        action='store_true',
        help='Disable automatic system configuration detection'
    )
    
    # Backend configuration
    backend_group = parser.add_argument_group('Backend Configuration')
    backend_group.add_argument(
        '--scheduler',
        choices=['slurm', 'pbs', 'local'],
        default='slurm',
        help='Job scheduler (default: slurm, or from system config)'
    )
    backend_group.add_argument(
        '--launcher',
        type=str,
        help='MPI launcher (auto-detected or from system config if not specified)'
    )
    backend_group.add_argument(
        '--module-system',
        choices=['lmod', 'tmod', 'tmod4', 'nomod'],
        default='lmod',
        help='Module system (default: lmod, or from system config)'
    )
    
    # Build configuration
    build_group = parser.add_argument_group('Build Configuration')
    build_group.add_argument(
        '--build-system',
        choices=['cmake', 'make', 'autotools'],
        help='Build system (auto-detected if not specified)'
    )
    build_group.add_argument(
        '--modules',
        type=str,
        help='Comma-separated list of modules to load (auto-detected if not specified)'
    )
    build_group.add_argument(
        '--no-build',
        action='store_true',
        help='Skip building, use existing executable'
    )
    
    # Input file configuration
    input_group = parser.add_argument_group('Input File Configuration')
    input_group.add_argument(
        '--input-file',
        type=str,
        help='Path to input file or directory (e.g., inputfiles/os-stdin)'
    )
    input_group.add_argument(
        '--input-name',
        type=str,
        help='Name of input file within input directory (e.g., os-stdin)'
    )
    input_group.add_argument(
        '--no-auto-input',
        action='store_true',
        help='Disable automatic input file detection'
    )
    
    # Behavior configuration
    behavior_group = parser.add_argument_group('Behavior Configuration')
    behavior_group.add_argument(
        '--no-submit',
        action='store_true',
        help='Generate job scripts but do not submit them'
    )
    behavior_group.add_argument(
        '--no-reports',
        action='store_true',
        help='Skip report generation'
    )
    behavior_group.add_argument(
        '--cleanup',
        action='store_true',
        help='Clean up cloned repositories after test'
    )
    behavior_group.add_argument(
        '--test-name',
        help='Custom test name (auto-generated if not specified)'
    )
    
    # Output configuration
    output_group = parser.add_argument_group('Output Configuration')
    output_group.add_argument(
        '--output-dir',
        type=Path,
        default=Path('output'),
        help='Output directory for results (default: output)'
    )
    output_group.add_argument(
        '--workspace-dir',
        type=Path,
        default=Path('workspace'),
        help='Workspace directory for Git clones (default: workspace)'
    )
    
    # Logging configuration
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug logging'
    )
    parser.add_argument(
        '--log-file',
        type=Path,
        help='Log to file instead of console'
    )
    
    return parser


def parse_tuple_arg(arg: str) -> tuple:
    """Parse comma-separated string into tuple of floats."""
    try:
        return tuple(float(x.strip()) for x in arg.split(','))
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid tuple format: {arg}. Expected format: x,y,z")


def load_system_configuration(args, logger) -> Optional[SystemConfigLoader]:
    """Load system configuration if available."""
    if args.no_system_config:
        logger.info("System configuration disabled by --no-system-config")
        return None
    
    # Try to load from specified file
    if args.system_config:
        if not args.system_config.exists():
            logger.warning(f"System config file not found: {args.system_config}")
            return None
        
        logger.info(f"Loading system configuration from {args.system_config}")
        loader = SystemConfigLoader(args.system_config)
        
        if loader.site_config:
            logger.info(f"✓ System configuration loaded")
            if loader.current_system:
                logger.info(f"✓ Detected system: {loader.current_system.name}")
            return loader
        else:
            logger.warning("Failed to load system configuration")
            return None
    
    # Try auto-detection
    logger.info("Attempting to auto-detect system configuration...")
    loader = auto_load_system_config()
    
    if loader and loader.current_system:
        logger.info(f"✓ Auto-detected system: {loader.current_system.name}")
        return loader
    
    logger.info("No system configuration found (continuing with manual settings)")
    return None


def apply_system_config_overrides(args, loader: Optional[SystemConfigLoader], logger):
    """Apply system configuration overrides to arguments."""
    if not loader:
        return args
    
    # Get partition name
    partition_name = args.partition_name
    if not partition_name and loader.current_system:
        # Use first non-login partition as default
        if loader.current_system.partitions:
            # Prefer compute partitions over login nodes
            for partition in loader.current_system.partitions:
                if partition.name != 'login' and partition.scheduler != 'local':
                    partition_name = partition.name
                    logger.info(f"Using default partition: {partition_name}")
                    break
            
            # Fallback to first partition if no compute partition found
            if not partition_name:
                partition_name = loader.current_system.partitions[0].name
                logger.info(f"Using default partition: {partition_name}")
    
    if not partition_name:
        logger.warning("No partition specified and no default available")
        return args
    
    # Get partition info
    partition = loader.get_partition(partition_name)
    if not partition:
        logger.warning(f"Partition '{partition_name}' not found in system config")
        return args
    
    logger.info(f"Applying system config from partition: {partition_name}")
    
    # Override hardware parameters if not manually specified
    if partition.processor:
        num_cpus = partition.processor.get('num_cpus', args.procs_per_node)
        if args.procs_per_node == 128:  # Default value
            args.procs_per_node = num_cpus
            logger.info(f"  procs_per_node: {num_cpus} (from system config)")
    
    if partition.devices and args.hardware == 'gpu':
        num_gpus = partition.devices[0].get('num_devices', 0)
        if args.gpus_per_node == 0:  # Default value
            args.gpus_per_node = num_gpus
            logger.info(f"  gpus_per_node: {num_gpus} (from system config)")
    
    # Override scheduler if not manually specified
    if args.scheduler == 'slurm':  # Default value
        args.scheduler = partition.scheduler
        logger.info(f"  scheduler: {partition.scheduler} (from system config)")
    
    # Override launcher if not manually specified and not using custom
    if not args.launcher and not args.use_custom_launcher:
        args.launcher = partition.launcher
        logger.info(f"  launcher: {partition.launcher} (from system config)")
    elif args.use_custom_launcher:
        args.launcher = args.use_custom_launcher
        logger.info(f"  launcher: {args.use_custom_launcher} (custom from system config)")
    
    # Override module system
    if loader.current_system and args.module_system == 'lmod':  # Default value
        args.module_system = loader.current_system.modules_system
        logger.info(f"  module_system: {loader.current_system.modules_system} (from system config)")
    
    # Override partition and account from access options
    if partition.access:
        for access in partition.access:
            if '--partition=' in access and args.partition == 'X_usr_prod':  # Default
                args.partition = access.split('=')[1]
                logger.info(f"  partition: {args.partition} (from system config)")
            elif '--account=' in access and args.account == 'cin_X':  # Default
                args.account = access.split('=')[1]
                logger.info(f"  account: {args.account} (from system config)")
    
    return args


def main():
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.debug else (logging.INFO if args.verbose else logging.WARNING)
    setup_logging(level=log_level, log_file=args.log_file)
    
    logger = logging.getLogger(__name__)
    logger.info("HPC Auto - Automated HPC Scaling Framework")
    logger.info("=" * 60)
    
    # Load system configuration
    system_loader = load_system_configuration(args, logger)
    
    # Apply system config overrides
    args = apply_system_config_overrides(args, system_loader, logger)
    
    logger.info("")  # Blank line for readability
    
    # Parse tuple arguments
    initial_procs = parse_tuple_arg(args.initial_procs)
    initial_domain = parse_tuple_arg(args.initial_domain) if args.initial_domain else None
    initial_cells = parse_tuple_arg(args.initial_cells) if args.initial_cells else None
    
    # Parse modules
    modules = [m.strip() for m in args.modules.split(',')] if args.modules else None
    
    # Create configuration
    config = OrchestratorConfig(
        source=args.source,
        scaling_type=args.scaling,
        hardware_type=args.hardware,
        max_nodes=args.nodes,
        initial_procs=initial_procs,
        initial_domain=initial_domain,
        initial_cells=initial_cells,
        procs_per_node=args.procs_per_node,
        gpus_per_node=args.gpus_per_node,
        time_limit=args.time_limit,
        partition=args.partition,
        account=args.account,
        scheduler=args.scheduler,
        launcher=args.launcher,
        module_system=args.module_system,
        build_system=args.build_system,
        modules=modules,
        input_file=args.input_file,
        input_file_name=args.input_name,
        auto_detect_input=not args.no_auto_input,
        auto_submit_jobs=not args.no_submit,
        cleanup_after_build=args.cleanup,
        generate_reports=not args.no_reports,
        output_dir=args.output_dir,
        workspace_dir=args.workspace_dir,
        test_name=args.test_name
    )
    
    # Create and run orchestrator
    orchestrator = HPCOrchestrator(config)
    success = orchestrator.run()
    
    # Exit with appropriate code
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
