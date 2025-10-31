#!/usr/bin/env python3
"""
Example: Using System Configuration for HPC Scaling Tests

This example demonstrates how to:
1. Load a system configuration file
2. Query partition information
3. Validate configurations
4. Use custom launchers
5. Generate scaling job configurations
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.system_loader import SystemConfigLoader
from engine.scaling import ScalingEngine
from core.config import ScalingConfig
from core.types import ScalingType
from core.factory import BackendFactory


def main():
    print("=" * 80)
    print("System Configuration Example")
    print("=" * 80)
    
    # Step 1: Load system configuration
    print("\n[1] Loading system configuration...")
    config_file = Path(__file__).parent.parent / "leonardo_system.py"
    
    if not config_file.exists():
        print(f"Error: Configuration file not found: {config_file}")
        print("Please ensure leonardo_system.py exists in the project root.")
        return 1
    
    loader = SystemConfigLoader(config_file)
    
    if not loader.site_config:
        print("Error: Failed to load system configuration")
        return 1
    
    print(f"✓ Loaded configuration from {config_file}")
    
    # Manually set current system if hostname doesn't match (for demo purposes)
    if not loader.current_system and loader.site_config.systems:
        print("  Note: Hostname doesn't match Leonardo patterns (expected for demo)")
        print("  Manually selecting 'leonardo' system for demonstration...")
        loader.current_system = loader.site_config.systems[0]
    
    # Step 2: Display system information
    print("\n[2] System Information:")
    loader.print_system_info()
    
    # Step 3: Select and validate partition
    print("\n[3] Partition Selection and Validation:")
    partition_name = "booster"
    partition = loader.get_partition(partition_name)
    
    if not partition:
        print(f"✗ Partition '{partition_name}' not found")
        return 1
    
    print(f"✓ Selected partition: {partition_name}")
    print(f"  Scheduler: {partition.scheduler}")
    print(f"  Launcher: {partition.launcher}")
    
    if partition.processor:
        print(f"  CPUs per node: {partition.processor.get('num_cpus', 'N/A')}")
    
    if partition.devices:
        gpu = partition.devices[0]
        print(f"  GPUs per node: {gpu.get('num_devices', 'N/A')}")
        print(f"  GPU Model: {gpu.get('model', 'N/A')}")
    
    # Step 4: Validate configuration
    print("\n[4] Configuration Validation:")
    
    # Import validation function from leonardo_system
    import leonardo_system
    
    test_nodes = 4
    test_procs_per_node = 32
    
    is_valid, message = leonardo_system.validate_scaling_config(
        partition_name=partition_name,
        num_nodes=test_nodes,
        procs_per_node=test_procs_per_node
    )
    
    if is_valid:
        print(f"✓ {message}")
    else:
        print(f"✗ {message}")
        return 1
    
    # Step 5: Create resource and backend configurations
    print("\n[5] Creating Resource Configuration:")
    
    resource_config = loader.create_resource_config(
        partition_name=partition_name,
        max_nodes=16,
        time_limit="01:00:00",
        account="MyProject"
    )
    
    print(f"  Max nodes: {resource_config.max_nodes}")
    print(f"  Procs per node: {resource_config.procs_per_node}")
    print(f"  GPUs per node: {resource_config.gpus_per_node}")
    print(f"  Time limit: {resource_config.time_limit}")
    print(f"  Partition: {resource_config.partition}")
    
    backend_config = loader.create_backend_config(
        partition_name=partition_name
    )
    
    print(f"\n  Backend Configuration:")
    print(f"  Scheduler: {backend_config.scheduler}")
    print(f"  Launcher: {backend_config.launcher}")
    print(f"  Module System: {backend_config.module_system}")
    
    # Step 6: Test custom launcher
    print("\n[6] Testing Custom Launchers:")
    
    from core.registry import list_launchers
    
    registered_launchers = list_launchers()
    print(f"  Registered custom launchers: {registered_launchers}")
    
    if 'mpirun-mapby' in registered_launchers:
        print("\n  Testing 'mpirun-mapby' launcher:")
        launcher = BackendFactory.create_launcher(
            'mpirun-mapby',
            options={'pe': 8}
        )
        
        # Create a test job config
        from core.config import JobConfig
        test_job = JobConfig(
            job_id="test",
            num_nodes=2,
            num_procs=64,
            procs_decomposition=(4, 4, 4)
        )
        
        launch_cmd = launcher.generate_launch_command(
            test_job,
            ['./my_app', '--input', 'data.in'],
            resource_config
        )
        
        print(f"  Generated command: {' '.join(launch_cmd)}")
    
    # Step 7: Generate scaling job configurations
    print("\n[7] Generating Scaling Job Configurations:")
    
    # Strong scaling configuration
    scaling_config = ScalingConfig(
        scaling_type=ScalingType.STRONG,
        max_nodes=16,
        initial_procs=(2, 2, 2),
        initial_domain=(10.0, 10.0, 10.0),
        initial_cells=(256, 256, 256)
    )
    
    print(f"  Scaling type: {scaling_config.scaling_type.value}")
    print(f"  Max nodes: {scaling_config.max_nodes}")
    print(f"  Initial decomposition: {scaling_config.initial_procs}")
    
    # Create scaling engine
    engine = ScalingEngine(scaling_config, resource_config)
    job_configs = engine.generate_job_configs()
    
    print(f"\n  Generated {len(job_configs)} job configurations:")
    print(f"  {'Job ID':<15} {'Nodes':<8} {'Procs':<8} {'Decomposition':<20} {'Domain Size'}")
    print("  " + "-" * 78)
    
    for job in job_configs:
        decomp_str = f"({job.procs_decomposition[0]}, {job.procs_decomposition[1]}, {job.procs_decomposition[2]})"
        domain_str = f"({job.domain_size[0]:.1f}, {job.domain_size[1]:.1f}, {job.domain_size[2]:.1f})" if job.domain_size else "N/A"
        
        print(f"  {job.job_id:<15} {job.num_nodes:<8} {job.num_procs:<8} {decomp_str:<20} {domain_str}")
    
    # Step 8: Display environment information
    print("\n[8] Available Environments:")
    
    if loader.site_config:
        for env in loader.site_config.environments:
            print(f"\n  {env.name}:")
            print(f"    Compilers: {env.cc}, {env.cxx}, {env.ftn}")
            if env.nvcc:
                print(f"    CUDA Compiler: {env.nvcc}")
            print(f"    Features: {', '.join(env.features)}")
            if env.modules:
                print(f"    Modules: {', '.join(env.modules)}")
    
    # Step 9: Summary
    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"""
✓ Successfully loaded system configuration
✓ Validated partition configuration
✓ Created resource and backend configurations
✓ Tested custom launcher registration
✓ Generated {len(job_configs)} scaling job configurations

The system is ready to run scaling tests!

Next steps:
1. Review the generated job configurations
2. Customize launcher options if needed
3. Run the scaling tests with: python hpc_auto.py <path> --nodes {resource_config.max_nodes}
    """)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
