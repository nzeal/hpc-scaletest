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
from utils.config_parser import load_yaml_config
from utils.system_info import auto_configure_resources, display_partition_info, get_partition_info


def create_parser():
    """Create argument parser."""
    parser = argparse.ArgumentParser(
        description="Automated HPC Scaling Framework - Simple, intelligent scaling tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Display system information
  python hpc_auto.py --system-info --system Leonardo --partition booster
  python hpc_auto.py --system-info --system Leonardo --partition dcgp

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
    
    # Required arguments (now optional if using --config)
    parser.add_argument(
        'source',
        nargs='?',  # Make optional
        help='Path to source code or Git repository URL (or use --config for YAML configuration)'
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
        '--system',
        type=str,
        help='System name (e.g., Leonardo) - requires system config file'
    )
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
        '--system-info',
        action='store_true',
        help='Display system information and exit (use with --system and --partition)'
    )
    behavior_group.add_argument(
        '--check-partition',
        action='store_true',
        help='Check partition hardware specs and exit (use with --partition)'
    )
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
        '--config',
        type=Path,
        help='Path to YAML configuration file (e.g., run.yaml)'
    )
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


def display_system_info_from_config(system_name: str, partition_name: str, system_config_path: Optional[Path] = None):
    """Display system information from system configuration file."""
    logger = logging.getLogger(__name__)
    
    # Try to find system config file
    if system_config_path and system_config_path.exists():
        config_file = system_config_path
    else:
        # Try standard location: <system_name>_system.py
        config_file = Path(f"{system_name.lower()}_system.py")
        if not config_file.exists():
            print(f"\n❌ Error: System configuration file not found: {config_file}")
            print(f"Please provide the path with --system-config or ensure {config_file.name} exists.")
            return False
    
    try:
        # Load system configuration
        loader = SystemConfigLoader(config_file)
        if not loader.site_config:
            print(f"\n❌ Error: Failed to load system configuration from {config_file}")
            return False
        
        # Find the system
        system = loader.site_config.get_system(system_name.lower())
        if not system:
            print(f"\n❌ Error: System '{system_name}' not found in configuration")
            available_systems = [s.name for s in loader.site_config.systems]
            if available_systems:
                print(f"Available systems: {', '.join(available_systems)}")
            return False
        
        # Find the partition
        partition = None
        for p in system.partitions:
            if p.name.lower() == partition_name.lower():
                partition = p
                break
        
        if not partition:
            print(f"\n❌ Error: Partition '{partition_name}' not found in system '{system_name}'")
            available_partitions = [p.name for p in system.partitions]
            if available_partitions:
                print(f"Available partitions: {', '.join(available_partitions)}")
            return False
        
        # Get partition information
        processor = partition.processor or {}
        devices = partition.devices or []
        
        # Display system information
        print("\n" + "="*80)
        print(f"SYSTEM INFORMATION: {system_name.upper()} - {partition_name.upper()} Partition")
        print("="*80)
        print(f"\nPartition Description: {partition.descr}")
        print(f"Scheduler:             {partition.scheduler}")
        print(f"Launcher:              {partition.launcher}")
        
        # CPU Information
        print("\n" + "-"*80)
        print("CPU CONFIGURATION")
        print("-"*80)
        cores_per_node = processor.get('num_cpus', 0)
        sockets = processor.get('num_sockets', 1)
        cores_per_socket = processor.get('num_cpus_per_socket', cores_per_node // sockets if sockets > 0 else cores_per_node)
        arch = processor.get('arch', 'Unknown')
        
        print(f"CPU Architecture:      {arch}")
        print(f"Cores per Node:        {cores_per_node}")
        print(f"Cores per Socket:      {cores_per_socket}")
        print(f"Sockets per Node:      {sockets}")
        
        # GPU Information
        if devices:
            print("\n" + "-"*80)
            print("GPU CONFIGURATION")
            print("-"*80)
            device = devices[0]  # Assume homogeneous GPUs
            gpu_count = device.get('num_devices', 0)
            gpu_model = device.get('model', 'Unknown')
            gpu_arch = device.get('arch', 'Unknown')
            
            print(f"GPU Model:             {gpu_model}")
            print(f"GPU Architecture:      {gpu_arch}")
            print(f"GPUs per Node:         {gpu_count}")
            
            # Try to get GPU memory from known models
            gpu_memory = {
                'A100': '40/80 GB',
                'V100': '16/32 GB',
                'H100': '80 GB',
                'A40': '48 GB'
            }.get(gpu_model, 'Check system specs')
            print(f"GPU Memory:            {gpu_memory}")
        else:
            print("\n" + "-"*80)
            print("GPU CONFIGURATION")
            print("-"*80)
            print("GPUs per Node:         0 (CPU-only partition)")
        
        # Try to get total nodes and memory from Slurm if available
        print("\n" + "-"*80)
        print("CLUSTER RESOURCES")
        print("-"*80)
        
        if partition.scheduler == 'slurm':
            try:
                # Extract partition name from access options
                slurm_partition = None
                for access_opt in partition.access:
                    if '--partition=' in access_opt:
                        slurm_partition = access_opt.split('=')[1]
                        break
                
                if slurm_partition:
                    print(f"Querying Slurm partition: {slurm_partition}\n")
                    # Use the existing get_partition_info function
                    partition_info = get_partition_info(slurm_partition)
                    
                    print(f"Total Nodes:           {partition_info['total_nodes']}")
                    print(f"Total Memory per Node: {partition_info['memory_per_node_gb']:.1f} GB")
                    
                    # Calculate totals
                    total_cores = partition_info['total_nodes'] * partition_info['cores_per_node']
                    total_memory = partition_info['total_nodes'] * partition_info['memory_per_node_gb']
                    
                    print(f"\nTotal CPU Cores:       {total_cores:,}")
                    print(f"Total Cluster Memory:  {total_memory:,.1f} GB ({total_memory/1024:.1f} TB)")
                    
                    if devices:
                        total_gpus = partition_info['total_nodes'] * partition_info['gpus_per_node']
                        print(f"Total GPUs:            {total_gpus:,}")
                else:
                    print("Note: Total nodes and memory require Slurm access")
                    print("      Run this command on the HPC system for complete info")
            except Exception as e:
                print(f"Note: Could not query Slurm (run on HPC system for full info)")
                logger.debug(f"Slurm query error: {e}")
        else:
            print("Note: Total nodes/memory info requires Slurm scheduler")
        
        # Additional information
        if partition.extras:
            print("\n" + "-"*80)
            print("ADDITIONAL INFORMATION")
            print("-"*80)
            for key, value in partition.extras.items():
                key_display = key.replace('_', ' ').title()
                print(f"{key_display:22s} {value}")
        
        print("="*80 + "\n")
        return True
        
    except Exception as e:
        print(f"\n❌ Error loading system configuration: {e}")
        import traceback
        traceback.print_exc()
        return False


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
    
    # Handle --system-info mode
    if args.system_info:
        if not args.system:
            logger.error("Please specify a system name with --system")
            logger.error("Example: python hpc_auto.py --system-info --system Leonardo --partition booster")
            sys.exit(1)
        
        if not args.partition and not args.partition_name:
            logger.error("Please specify a partition with --partition or --partition-name")
            logger.error("Example: python hpc_auto.py --system-info --system Leonardo --partition booster")
            sys.exit(1)
        
        partition = args.partition_name if args.partition_name else args.partition
        if partition == 'X_usr_prod':  # Default value
            logger.error("Please specify a valid partition name")
            logger.error("Example: python hpc_auto.py --system-info --system Leonardo --partition booster")
            sys.exit(1)
        
        success = display_system_info_from_config(args.system, partition, args.system_config)
        sys.exit(0 if success else 1)
    
    # Handle --check-partition mode
    if args.check_partition:
        if not args.partition or args.partition == 'X_usr_prod':
            logger.error("Please specify a partition name with --partition")
            logger.error("Example: python hpc_auto.py --check-partition --partition dcgp")
            sys.exit(1)
        
        logger.info(f"Checking partition: {args.partition}")
        display_partition_info(args.partition)
        sys.exit(0)
    
    # Check if using YAML configuration
    if args.config:
        logger.info(f"Loading configuration from {args.config}")
        try:
            # First, check if YAML has verbose flag to set up logging early
            import yaml
            with open(args.config, 'r') as f:
                yaml_data = yaml.safe_load(f)
            
            # Re-setup logging if YAML specifies verbose
            if yaml_data and yaml_data.get('verbose', False):
                if not args.verbose and not args.debug:
                    # YAML has verbose=true, update logging
                    setup_logging(level=logging.INFO, log_file=args.log_file)
                    logger.setLevel(logging.INFO)
            
            # Now load full YAML config
            yaml_config = load_yaml_config(args.config)
            
            # Merge with command-line arguments (CLI takes precedence)
            config_dict = yaml_config.copy()
            
            # Override with CLI arguments if provided
            if args.source:
                config_dict['source'] = args.source
            
            # Check if source is specified either way
            if 'source' not in config_dict:
                logger.error("No source specified in YAML config or command line")
                sys.exit(1)
            
            # Apply other CLI overrides
            if args.nodes != 4:  # Not default
                config_dict['max_nodes'] = args.nodes
            if args.scaling != 'strong':  # Not default
                config_dict['scaling_type'] = args.scaling
            if args.partition != 'X_usr_prod':  # Not default
                config_dict['partition'] = args.partition
            if args.account != 'cin_X':  # Not default
                config_dict['account'] = args.account
            
            # Auto-detect system resources if not explicitly configured
            logger.info("\nDetecting system resources...")
            partition_for_detection = config_dict.get('partition', args.partition)
            if partition_for_detection == 'X_usr_prod':  # Default placeholder
                partition_for_detection = None
            
            auto_config = auto_configure_resources(
                max_nodes=config_dict.get('max_nodes'),
                partition=partition_for_detection
            )
            
            # Apply detected values if not explicitly set in YAML
            if 'procs_per_node' not in yaml_config:
                config_dict['procs_per_node'] = auto_config['procs_per_node']
                logger.info(f"  Auto-detected procs_per_node: {auto_config['procs_per_node']}")
            
            if 'gpus_per_node' not in yaml_config:
                config_dict['gpus_per_node'] = auto_config['gpus_per_node']
                if auto_config['gpus_per_node'] > 0:
                    logger.info(f"  Auto-detected gpus_per_node: {auto_config['gpus_per_node']}")
            
            if 'memory_per_node_gb' not in yaml_config:
                config_dict['memory_per_node_gb'] = auto_config['memory_per_node_gb']
                logger.info(f"  Auto-detected memory: {auto_config['memory_per_node_gb']:.1f} GB/node")
            
            if 'scheduler' not in yaml_config:
                config_dict['scheduler'] = auto_config['scheduler']
                logger.info(f"  Auto-detected scheduler: {auto_config['scheduler']}")
            
            # Use resolved partition name if auto-detected
            if auto_config.get('partition') and partition_for_detection:
                config_dict['partition'] = auto_config['partition']
                if auto_config['partition'] != partition_for_detection:
                    logger.info(f"  Resolved partition '{partition_for_detection}' to '{auto_config['partition']}'")
            
            # Create configuration from merged settings
            config = OrchestratorConfig(**config_dict)
            
            logger.info("Configuration loaded from YAML")
            
        except Exception as e:
            logger.error(f"Failed to load YAML configuration: {e}")
            sys.exit(1)
    else:
        # Original CLI-based configuration
        if not args.source:
            logger.error("Error: source argument is required when not using --config")
            parser.print_help()
            sys.exit(1)
        
        # Load system configuration
        system_loader = load_system_configuration(args, logger)
        
        # Apply system config overrides
        args = apply_system_config_overrides(args, system_loader, logger)
        
        # Auto-detect system resources if not explicitly set
        logger.info("\nDetecting system resources...")
        partition_for_detection = args.partition if args.partition != 'X_usr_prod' else None
        
        auto_config = auto_configure_resources(
            max_nodes=args.nodes,
            partition=partition_for_detection
        )
        
        # Apply detected values if using defaults (not explicitly set by user)
        if args.procs_per_node == 128:  # Default value
            args.procs_per_node = auto_config['procs_per_node']
            logger.info(f"  Auto-detected procs_per_node: {auto_config['procs_per_node']}")
        
        if args.gpus_per_node == 0 and args.hardware == 'gpu':  # Default for GPU
            args.gpus_per_node = auto_config['gpus_per_node']
            if auto_config['gpus_per_node'] > 0:
                logger.info(f"  Auto-detected gpus_per_node: {auto_config['gpus_per_node']}")
        
        # Always detect memory
        memory_gb = auto_config['memory_per_node_gb']
        logger.info(f"  Auto-detected memory: {memory_gb:.1f} GB/node")
        
        if args.scheduler == 'slurm':  # Default value
            args.scheduler = auto_config['scheduler']
            logger.info(f"  Auto-detected scheduler: {auto_config['scheduler']}")
        
        logger.info("")  # Blank line for readability
    
        # Parse tuple arguments
        initial_procs = parse_tuple_arg(args.initial_procs)
        initial_domain = parse_tuple_arg(args.initial_domain) if args.initial_domain else None
        initial_cells = parse_tuple_arg(args.initial_cells) if args.initial_cells else None
        
        # CRITICAL FIX: If initial_procs is provided, calculate procs_per_node from it
        # This ensures weak scaling matches the MPI decomposition specified by user
        if initial_procs:
            calculated_procs_per_node = initial_procs[0] * initial_procs[1] * initial_procs[2]
            logger.info(f"  Calculated procs_per_node from initial_procs: {initial_procs} = {calculated_procs_per_node}")
            args.procs_per_node = calculated_procs_per_node  # Override with calculated value
        elif args.procs_per_node == 128:  # Default value, not explicitly set
            args.procs_per_node = auto_config['procs_per_node']
            logger.info(f"  Auto-detected procs_per_node: {auto_config['procs_per_node']}")
        
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
            memory_per_node_gb=memory_gb,
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
