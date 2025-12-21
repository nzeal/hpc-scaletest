#!/usr/bin/env python3
"""
Diagnostic tool for HPC-ScaleTest GPU job configuration.
Previews generated job scripts and validates configuration before submission.

Usage:
    python diagnose_gpu_config.py --config runGpuWeakScaling.yaml --nodes 1,2,4
"""

import argparse
import yaml
import sys
from pathlib import Path

def load_config(config_file):
    """Load YAML configuration file."""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return config

def analyze_gpu_config(config):
    """Analyze GPU configuration and check for issues."""
    print("=" * 80)
    print("GPU CONFIGURATION ANALYSIS")
    print("=" * 80)
    
    hardware_type = config.get('hardware_type', 'unknown')
    gpus_per_node = config.get('gpus_per_node', 0)
    procs_per_node = config.get('procs_per_node', 0)
    initial_procs = config.get('initial_procs', [1, 1, 1])
    
    print(f"\n✓ Hardware Type: {hardware_type}")
    print(f"✓ GPUs per Node: {gpus_per_node}")
    print(f"✓ Procs per Node: {procs_per_node}")
    print(f"✓ Initial Procs: {initial_procs[0]} × {initial_procs[1]} × {initial_procs[2]} = "
          f"{initial_procs[0] * initial_procs[1] * initial_procs[2]}")
    
    # Validate configuration
    issues = []
    warnings = []
    
    if hardware_type == 'gpu':
        total_initial_procs = initial_procs[0] * initial_procs[1] * initial_procs[2]
        
        # Check 1: initial_procs should match gpus_per_node for weak scaling
        if total_initial_procs != gpus_per_node:
            issues.append(
                f"⚠ MISMATCH: initial_procs ({total_initial_procs}) != gpus_per_node ({gpus_per_node})\n"
                f"   For weak scaling, initial_procs should equal gpus_per_node (1 MPI task per GPU)"
            )
        
        # Check 2: procs_per_node should equal gpus_per_node for GPU jobs
        if procs_per_node != gpus_per_node:
            if procs_per_node > gpus_per_node:
                # Might be specifying CPU cores instead of MPI tasks
                warnings.append(
                    f"⚠ INFO: procs_per_node ({procs_per_node}) > gpus_per_node ({gpus_per_node})\n"
                    f"   This will be interpreted as: {procs_per_node} CPU cores per node\n"
                    f"   MPI tasks will be set to {gpus_per_node} (1 per GPU)\n"
                    f"   Each task gets ~{procs_per_node // gpus_per_node} CPU cores"
                )
            else:
                issues.append(
                    f"⚠ ERROR: procs_per_node ({procs_per_node}) < gpus_per_node ({gpus_per_node})\n"
                    f"   Not enough MPI tasks for available GPUs!"
                )
    
    # Print issues
    if issues:
        print("\n" + "!" * 80)
        print("CONFIGURATION ISSUES DETECTED:")
        print("!" * 80)
        for issue in issues:
            print(f"\n{issue}")
    
    if warnings:
        print("\n" + "-" * 80)
        print("CONFIGURATION WARNINGS:")
        print("-" * 80)
        for warning in warnings:
            print(f"\n{warning}")
    
    if not issues and not warnings:
        print("\n✓ Configuration looks good!")
    
    return len(issues) == 0

def preview_job_script(config, num_nodes):
    """Preview what the job script would look like for given node count."""
    print("\n" + "=" * 80)
    print(f"JOB SCRIPT PREVIEW FOR {num_nodes} NODE(S)")
    print("=" * 80)
    
    hardware_type = config.get('hardware_type', 'cpu')
    gpus_per_node = config.get('gpus_per_node', 0)
    procs_per_node = config.get('procs_per_node', 1)
    initial_procs = config.get('initial_procs', [1, 1, 1])
    partition = config.get('partition', 'unknown')
    qos = config.get('qos', 'unknown')
    account = config.get('account', 'unknown')
    time_limit = config.get('time_limit', '01:00:00')
    
    # Calculate MPI decomposition for this node count
    total_initial_procs = initial_procs[0] * initial_procs[1] * initial_procs[2]
    
    if hardware_type == 'gpu':
        # GPU job
        tasks_per_node = gpus_per_node
        total_tasks = num_nodes * tasks_per_node
        
        # Calculate CPUs per task
        if procs_per_node > gpus_per_node:
            # procs_per_node is CPU cores
            cpus_per_task = procs_per_node // gpus_per_node
        else:
            # Assume 1 CPU per task as fallback
            cpus_per_task = 1
        
        # MPI decomposition (for weak scaling, scales in 2D)
        scaling_dims = config.get('scaling_dimensions', 2)
        mpi_x = initial_procs[0]
        mpi_y = initial_procs[1]
        mpi_z = initial_procs[2]
        
        # Scale up based on node count
        factor = num_nodes
        if scaling_dims == 2:
            # Scale in X and Y alternately
            import math
            sqrt_factor = int(math.sqrt(factor))
            if sqrt_factor * sqrt_factor == factor:
                mpi_x *= sqrt_factor
                mpi_y *= sqrt_factor
            else:
                # Scale X first, then Y
                while factor > 1:
                    if mpi_x <= mpi_y:
                        mpi_x *= 2
                    else:
                        mpi_y *= 2
                    factor //= 2
        
        print("\nSLURM Directives:")
        print(f"  #SBATCH --nodes={num_nodes}")
        print(f"  #SBATCH --partition={partition}")
        print(f"  #SBATCH --qos={qos}")
        print(f"  #SBATCH --ntasks={total_tasks}")
        print(f"  #SBATCH --ntasks-per-node={tasks_per_node}")
        print(f"  #SBATCH --cpus-per-task={cpus_per_task}")
        print(f"  #SBATCH --gres=gpu:{gpus_per_node}")
        print(f"  #SBATCH -A {account}")
        print(f"  #SBATCH --time={time_limit}")
        
        print("\nMPI Launch Command:")
        print(f"  mpirun -np {total_tasks} \\")
        print(f"    --map-by ppr:{tasks_per_node}:node \\")
        print(f"    --bind-to core \\")
        print(f"    --cpus-per-proc {cpus_per_task} \\")
        print(f"    ./bind_gpu $BINARY/executable input_file")
        
        print("\nResource Summary:")
        print(f"  • Nodes: {num_nodes}")
        print(f"  • GPUs per node: {gpus_per_node}")
        print(f"  • Total GPUs: {num_nodes * gpus_per_node}")
        print(f"  • MPI tasks per node: {tasks_per_node} (1 per GPU)")
        print(f"  • Total MPI tasks: {total_tasks}")
        print(f"  • CPUs per MPI task: {cpus_per_task}")
        print(f"  • MPI decomposition: {mpi_x} × {mpi_y} × {mpi_z} = {mpi_x * mpi_y * mpi_z}")
        
    else:
        # CPU job
        total_tasks = num_nodes * procs_per_node
        
        print("\nSLURM Directives:")
        print(f"  #SBATCH --nodes={num_nodes}")
        print(f"  #SBATCH --partition={partition}")
        print(f"  #SBATCH --qos={qos}")
        print(f"  #SBATCH --ntasks={total_tasks}")
        print(f"  #SBATCH --ntasks-per-node={procs_per_node}")
        print(f"  #SBATCH --cpus-per-task=1")
        print(f"  #SBATCH -A {account}")
        print(f"  #SBATCH --time={time_limit}")
        
        print("\nMPI Launch Command:")
        print(f"  mpirun -np {total_tasks} \\")
        print(f"    --map-by ppr:{procs_per_node}:node \\")
        print(f"    $BINARY/executable input_file")
        
        print("\nResource Summary:")
        print(f"  • Nodes: {num_nodes}")
        print(f"  • MPI tasks per node: {procs_per_node}")
        print(f"  • Total MPI tasks: {total_tasks}")

def check_qos_requirements(config, num_nodes):
    """Check if configuration meets QOS requirements."""
    print("\n" + "=" * 80)
    print("QOS REQUIREMENTS CHECK")
    print("=" * 80)
    
    partition = config.get('partition', '')
    qos = config.get('qos', '')
    gpus_per_node = config.get('gpus_per_node', 0)
    procs_per_node = config.get('procs_per_node', 1)
    
    print(f"\n✓ Partition: {partition}")
    print(f"✓ QOS: {qos}")
    
    # Known QOS requirements for Leonardo Booster
    if 'boost' in partition.lower() and 'boost' in qos.lower():
        print("\n📋 Leonardo Booster QOS Requirements:")
        print("  • boost_qos_bprod typically requires:")
        print("    - Minimum 2 nodes OR minimum 64 CPUs")
        print("    - GPU jobs must request GPUs with --gres=gpu:N")
        
        # Check if requirements are met
        if gpus_per_node > 0:
            # GPU job
            tasks_per_node = gpus_per_node
            cpus_per_task = procs_per_node // gpus_per_node if procs_per_node > gpus_per_node else 1
            total_cpus = num_nodes * tasks_per_node * cpus_per_task
        else:
            # CPU job
            total_cpus = num_nodes * procs_per_node
        
        print(f"\n✓ Your configuration:")
        print(f"  • Nodes: {num_nodes}")
        print(f"  • Estimated total CPUs: {total_cpus}")
        
        if num_nodes == 1 and total_cpus < 64:
            print("\n⚠ WARNING: Single-node job with < 64 CPUs may not meet QOS minimum!")
            print("  Recommended solutions:")
            print("    1. Use at least 2 nodes")
            print("    2. Or use a different QOS (e.g., boost_qos_dbg for debugging)")
            print("    3. Or increase cpus-per-task to reach 64 total CPUs")
        else:
            print("\n✓ Configuration likely meets QOS requirements")
    else:
        print("\n📋 QOS requirements unknown for this partition/qos combination")
        print("  Please consult your HPC system documentation")

def main():
    parser = argparse.ArgumentParser(
        description="Diagnose HPC-ScaleTest GPU configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze configuration and preview for 1, 2, and 4 nodes
  python diagnose_gpu_config.py --config runGpuWeakScaling.yaml --nodes 1,2,4
  
  # Quick check without preview
  python diagnose_gpu_config.py --config runGpuWeakScaling.yaml --check-only
        """
    )
    parser.add_argument('--config', required=True, help='Path to YAML configuration file')
    parser.add_argument('--nodes', help='Comma-separated list of node counts to preview (e.g., 1,2,4,8)')
    parser.add_argument('--check-only', action='store_true', help='Only validate configuration, no preview')
    
    args = parser.parse_args()
    
    # Load configuration
    try:
        config = load_config(args.config)
    except Exception as e:
        print(f"❌ Error loading config file: {e}")
        sys.exit(1)
    
    # Analyze configuration
    is_valid = analyze_gpu_config(config)
    
    if args.check_only:
        sys.exit(0 if is_valid else 1)
    
    # Preview job scripts for specified node counts
    if args.nodes:
        node_counts = [int(n.strip()) for n in args.nodes.split(',')]
    else:
        # Default preview for 1 node
        node_counts = [1]
    
    for num_nodes in node_counts:
        preview_job_script(config, num_nodes)
        check_qos_requirements(config, num_nodes)
    
    print("\n" + "=" * 80)
    print("DIAGNOSIS COMPLETE")
    print("=" * 80)
    
    if not is_valid:
        print("\n⚠ Configuration issues detected. Please review and fix before running.")
        sys.exit(1)
    else:
        print("\n✓ Configuration validated successfully!")
        print("  You can now run: python hpc_auto.py --config", args.config)
        sys.exit(0)

if __name__ == '__main__':
    main()
