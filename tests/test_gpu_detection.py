#!/usr/bin/env python3
"""
Simple GPU Detection Test - Python 3.6 Compatible

Tests partition detection with flexible naming.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.advanced_gpu_manager import AdvancedGPUManager
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')

if len(sys.argv) < 2:
    print("Usage: python3 test_gpu_detection.py <partition>")
    print("")
    print("Examples:")
    print("  python3 test_gpu_detection.py booster")
    print("  python3 test_gpu_detection.py boost_usr_prod")
    print("  python3 test_gpu_detection.py boost_qos_dbg")
    print("")
    print("The script will automatically find the correct partition name.")
    sys.exit(1)

partition = sys.argv[1]

print("="*70)
print("Testing GPU Detection with Partition: '{}'".format(partition))
print("="*70)
print("")

try:
    manager = AdvancedGPUManager()
    config = manager.detect_gpu_node_config(partition)
    
    if config:
        print("")
        print("="*70)
        print("SUCCESS! GPU Configuration Detected:")
        print("="*70)
        print("  Actual partition:  {}".format(config.partition))
        print("  GPUs per node:     {}".format(config.gpus_per_node))
        print("  GPU model:         {}".format(config.gpu_model))
        print("  CPUs per node:     {}".format(config.cpus_per_node))
        print("  Cores per GPU:     {}".format(config.cores_per_gpu))
        print("  Memory:            {} GB/node".format(config.memory_per_node_gb))
        print("")
        
        print("="*70)
        print("Example Job Configuration:")
        print("="*70)
        print("")
        print("SLURM Directives:")
        directives = manager.generate_slurm_directives(num_nodes=1)
        print("  #SBATCH --nodes=1")
        print("  #SBATCH --partition={}".format(directives['partition']))
        print("  #SBATCH --ntasks-per-node={}    # CPU allocation".format(directives['ntasks-per-node']))
        print("  #SBATCH --gres=gpu:{}             # GPU allocation".format(config.gpus_per_node))
        print("")
        
        print("Launch Commands:")
        print("")
        print("  Option A (mpirun):")
        mpirun = manager.generate_mpirun_command("$BINARY/app", "input.txt")
        print("    {}".format(mpirun))
        print("")
        print("  Option B (srun):")
        srun = manager.generate_srun_command("$BINARY/app", "input.txt")
        print("    {}".format(srun))
        print("")
        
        print("="*70)
        print("")
        
    else:
        print("")
        print("❌ Could not detect GPU configuration for partition '{}'".format(partition))
        print("")
        print("Troubleshooting:")
        print("  1. Check if partition exists: sinfo -o '%R' | grep {}".format(partition))
        print("  2. Check partition details: scontrol show partition {}".format(partition))
        print("  3. Make sure you're on a login node with SLURM access")
        print("")
        sys.exit(1)
        
except Exception as e:
    print("")
    print("❌ Error during detection: {}".format(e))
    print("")
    import traceback
    traceback.print_exc()
    sys.exit(1)
