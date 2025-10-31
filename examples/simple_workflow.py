#!/usr/bin/env python3
"""
Simple workflow example for HPC Auto.
Demonstrates the easiest way to run scaling tests.
"""

from pathlib import Path
from engine.orchestrator import create_simple_workflow
from utils.logging_config import setup_logging


def main():
    """Run a simple scaling test."""
    # Setup logging
    setup_logging()
    
    # Example 1: Strong scaling test on local code
    print("\n" + "="*60)
    print("Example 1: Strong Scaling Test - Local Code")
    print("="*60)
    
    orchestrator = create_simple_workflow(
        source="/path/to/your/hpc/code",  # Change this to your code path
        scaling_type="strong",
        max_nodes=8,
        hardware_type="cpu"
    )
    
    # Run the workflow
    success = orchestrator.run()
    print(f"\nWorkflow completed: {'SUCCESS' if success else 'FAILED'}")
    
    
    # Example 2: Weak scaling test from Git
    print("\n" + "="*60)
    print("Example 2: Weak Scaling Test - Git Repository")
    print("="*60)
    
    orchestrator = create_simple_workflow(
        source="https://github.com/user/hpc-app.git",  # Your Git URL
        scaling_type="weak",
        max_nodes=16,
        hardware_type="cpu",
        initial_domain=(10.0, 10.0, 10.0),  # Initial domain size for weak scaling
        initial_cells=(256, 256, 256)  # Initial cell count
    )
    
    success = orchestrator.run()
    print(f"\nWorkflow completed: {'SUCCESS' if success else 'FAILED'}")
    
    
    # Example 3: GPU-enabled test
    print("\n" + "="*60)
    print("Example 3: GPU Scaling Test")
    print("="*60)
    
    orchestrator = create_simple_workflow(
        source="/path/to/cuda/app",
        scaling_type="strong",
        max_nodes=4,
        hardware_type="gpu",
        procs_per_node=4,  # Fewer processes for GPU nodes
        time_limit="01:00:00"
    )
    
    success = orchestrator.run()
    print(f"\nWorkflow completed: {'SUCCESS' if success else 'FAILED'}")


if __name__ == '__main__':
    main()
