#!/usr/bin/env python3
"""
Test GPU Scaling Configuration

Verifies that GPU task layout is correctly calculated for:
- Leonardo: 4 GPUs per node, 112 cores per node
- Generic systems with different GPU counts
"""

import sys
sys.path.insert(0, '.')

from core.config import ResourceConfig

def test_gpu_task_configuration():
    """Test GPU task calculation for various configurations."""
    
    print("="*70)
    print("GPU TASK CONFIGURATION TESTS")
    print("="*70)
    print()
    
    # Test 1: Leonardo configuration (4 GPUs, 112 cores)
    print("Test 1: Leonardo (4 GPUs, 112 cores)")
    print("-"*70)
    
    config1 = ResourceConfig(
        max_nodes=4,
        procs_per_node=112,
        gpus_per_node=4
    )
    
    config1.configure_gpu_tasks(cores_per_node=112)
    
    assert config1.actual_mpi_tasks == 4, f"Expected 4 tasks, got {config1.actual_mpi_tasks}"
    assert config1.cores_per_task == 28, f"Expected 28 cores/task, got {config1.cores_per_task}"
    
    print(f"✅ MPI tasks per node: {config1.actual_mpi_tasks}")
    print(f"✅ CPU cores per task: {config1.cores_per_task}")
    print(f"✅ For 1 node: {config1.actual_mpi_tasks * 1} total MPI tasks")
    print(f"✅ For 2 nodes: {config1.actual_mpi_tasks * 2} total MPI tasks")
    print(f"✅ For 4 nodes: {config1.actual_mpi_tasks * 4} total MPI tasks")
    print()
    
    # Test 2: NVIDIA DGX A100 (8 GPUs, 128 cores)
    print("Test 2: NVIDIA DGX A100 (8 GPUs, 128 cores)")
    print("-"*70)
    
    config2 = ResourceConfig(
        max_nodes=2,
        procs_per_node=128,
        gpus_per_node=8
    )
    
    config2.configure_gpu_tasks(cores_per_node=128)
    
    assert config2.actual_mpi_tasks == 8, f"Expected 8 tasks, got {config2.actual_mpi_tasks}"
    assert config2.cores_per_task == 16, f"Expected 16 cores/task, got {config2.cores_per_task}"
    
    print(f"✅ MPI tasks per node: {config2.actual_mpi_tasks}")
    print(f"✅ CPU cores per task: {config2.cores_per_task}")
    print(f"✅ For 2 nodes: {config2.actual_mpi_tasks * 2} total MPI tasks")
    print()
    
    # Test 3: CPU-only (no GPUs)
    print("Test 3: CPU-only configuration (no GPUs)")
    print("-"*70)
    
    config3 = ResourceConfig(
        max_nodes=4,
        procs_per_node=128,
        gpus_per_node=0
    )
    
    config3.configure_gpu_tasks(cores_per_node=128)
    
    assert config3.actual_mpi_tasks is None, "CPU-only should not set actual_mpi_tasks"
    assert config3.cores_per_task is None, "CPU-only should not set cores_per_task"
    
    print(f"✅ CPU-only: No GPU task configuration")
    print(f"✅ Uses procs_per_node ({config3.procs_per_node}) for MPI tasks")
    print()
    
    # Test 4: Single GPU systems
    print("Test 4: Single GPU system (1 GPU, 32 cores)")
    print("-"*70)
    
    config4 = ResourceConfig(
        max_nodes=8,
        procs_per_node=32,
        gpus_per_node=1
    )
    
    config4.configure_gpu_tasks(cores_per_node=32)
    
    assert config4.actual_mpi_tasks == 1, f"Expected 1 task, got {config4.actual_mpi_tasks}"
    assert config4.cores_per_task == 32, f"Expected 32 cores/task, got {config4.cores_per_task}"
    
    print(f"✅ MPI tasks per node: {config4.actual_mpi_tasks}")
    print(f"✅ CPU cores per task: {config4.cores_per_task}")
    print(f"✅ For 8 nodes: {config4.actual_mpi_tasks * 8} total MPI tasks")
    print()
    
    print("="*70)
    print("✅ ALL GPU TASK CONFIGURATION TESTS PASSED!")
    print("="*70)


def test_mpi_command_structure():
    """Verify MPI command structure for GPU runs."""
    
    print()
    print("="*70)
    print("EXPECTED MPI COMMAND STRUCTURES (CORRECT SYNTAX)")
    print("="*70)
    print()
    
    # Leonardo example
    print("Leonardo (4 GPUs, 112 cores, 2 nodes):")
    print("-"*70)
    config = ResourceConfig(procs_per_node=112, gpus_per_node=4)
    config.configure_gpu_tasks(112)
    
    total_tasks = config.actual_mpi_tasks * 2  # 2 nodes
    
    # FIXED: Correct OpenMPI syntax with separate --map-by and --bind-to
    print(f"  mpirun -np {total_tasks} \\")
    print(f"    --map-by ppr:{config.actual_mpi_tasks}:node \\")
    print(f"    --bind-to core \\")
    print(f"    --cpus-per-proc {config.cores_per_task} \\")
    print(f"    --report-bindings \\")
    print(f"    ./app")
    print()
    print(f"  Explanation:")
    print(f"    -np {total_tasks}: {total_tasks} total MPI tasks (4 per node × 2 nodes)")
    print(f"    --map-by ppr:4:node: 4 processes per node (1 per GPU)")
    print(f"    --bind-to core: Bind processes to CPU cores")
    print(f"    --cpus-per-proc 28: 28 CPU cores per process")
    print()
    
    # DGX example
    print("NVIDIA DGX A100 (8 GPUs, 128 cores, 1 node):")
    print("-"*70)
    config2 = ResourceConfig(procs_per_node=128, gpus_per_node=8)
    config2.configure_gpu_tasks(128)
    
    total_tasks2 = config2.actual_mpi_tasks * 1  # 1 node
    
    print(f"  mpirun -np {total_tasks2} \\")
    print(f"    --map-by ppr:{config2.actual_mpi_tasks}:node \\")
    print(f"    --bind-to core \\")
    print(f"    --cpus-per-proc {config2.cores_per_task} \\")
    print(f"    --report-bindings \\")
    print(f"    ./app")
    print()
    print(f"  Explanation:")
    print(f"    -np {total_tasks2}: {total_tasks2} total MPI tasks (8 per node × 1 node)")
    print(f"    --map-by ppr:8:node: 8 processes per node (1 per GPU)")
    print(f"    --bind-to core: Bind processes to CPU cores")
    print(f"    --cpus-per-proc 16: 16 CPU cores per process")
    print()
    
    print("="*70)


if __name__ == "__main__":
    test_gpu_task_configuration()
    test_mpi_command_structure()
    print()
    print("✅ All GPU scaling tests PASSED!")
    print()
