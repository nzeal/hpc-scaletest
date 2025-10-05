#!/usr/bin/env python3
"""
Complete working example for HPC-ScaleTest.

This example demonstrates a realistic CFD benchmark scaling study
that can run both locally (for testing) and on a Slurm cluster.

To run this example:
1. Local testing:
   python scaletest.py run --test complete_example.py --backend local --max-nodes 1

2. On Slurm cluster:
   python scaletest.py run --test complete_example.py --max-nodes 64

3. Validation:
   python scaletest.py validate --test complete_example.py
"""

from pathlib import Path
from core.test_definition import Test


def create_cfd_strong_scaling_test():
    """
    Create a comprehensive CFD strong scaling test.
    
    This test simulates a computational fluid dynamics benchmark
    with a fixed problem size (512^3 cells) scaling from 1 to 128 nodes.
    
    Returns:
        Test: Configured test instance
    """
    
    # Create test instance
    test = Test(
        name="cfd_strong_scaling",
        input_file=Path("./inputs/cfd_512cubed.inp"),
        command=["./cfd_solver", "--solver", "gmres", "--tolerance", "1e-6"]
    )
    
    # Configure environment modules
    # These modules provide the compiler toolchain and MPI implementation
    test.set_modules([
        "gcc/11.2.0",           # GNU compiler
        "openmpi/4.1.1",        # MPI library
        "hdf5/1.12.1",          # For I/O
        "fftw/3.3.10"           # FFT library
    ])
    
    # Set environment variables
    # These tune the MPI and application behavior
    test.set_env({
        # OpenMP threading (disable for pure MPI)
        "OMP_NUM_THREADS": "1",
        
        # MPI tuning
        "MPI_TYPE_DEPTH": "16",
        "OMPI_MCA_btl": "^openib",  # Disable deprecated OpenIB
        
        # I/O optimization
        "HDF5_USE_FILE_LOCKING": "FALSE",
        "ROMIO_HINTS": "./romio_hints.txt",
        
        # Application-specific
        "CFD_OUTPUT_FREQ": "100",
        "CFD_CHECKPOINT_FREQ": "500"
    })
    
    # Configure backend components
    # This determines how jobs are submitted and executed
    test.set_backend(
        scheduler="slurm",       # Use Slurm workload manager
        launcher="srun",         # Use srun for parallel execution
        module_system="lmod",    # Use Lmod for modules
        build_system="cmake"     # Use CMake if building
    )
    
    # Configure computational resources
    # Adjust these based on your cluster configuration
    test.set_resources(
        max_nodes=128,                  # Scale up to 128 nodes
        procs_per_node=128,             # 128 MPI ranks per node
        gpus_per_node=0,                # CPU-only job
        memory_per_node="100GB",        # Memory per node
        time_limit="04:00:00",          # 4 hour time limit
        partition="compute",            # Compute partition
        account="cfd_project",          # Billing account
    )
    
    # Configure scaling parameters
    # This defines the problem decomposition and scaling strategy
    test.set_scaling(
        scaling_type="strong",          # Fixed problem size
        max_nodes=128,
        
        # Initial 3D processor decomposition: 4x4x4 = 64 processes
        # This will be doubled alternately in x and y directions
        initial_procs=(4, 4, 4),
        
        # Problem domain size (constant for strong scaling)
        # 100 x 100 x 100 physical domain
        initial_domain=(100.0, 100.0, 100.0),
        
        # Grid resolution (constant for strong scaling)
        # 512 x 512 x 512 = 134 million cells
        initial_cells=(512, 512, 512)
    )
    
    return test


def create_cfd_weak_scaling_test():
    """
    Create a CFD weak scaling test.
    
    This test maintains constant work per processor by growing
    the problem size proportionally with processor count.
    
    Returns:
        Test: Configured test instance
    """
    
    test = Test(
        name="cfd_weak_scaling",
        input_file=Path("./inputs/cfd_weak.inp"),
        command=["./cfd_solver", "--solver", "gmres"]
    )
    
    test.set_modules(["gcc/11.2.0", "openmpi/4.1.1", "hdf5/1.12.1"])
    test.set_env({"OMP_NUM_THREADS": "1", "HDF5_USE_FILE_LOCKING": "FALSE"})
    
    test.set_backend(
        scheduler="slurm",
        launcher="srun",
        module_system="lmod"
    )
    
    test.set_resources(
        max_nodes=256,
        procs_per_node=128,
        time_limit="08:00:00",
        partition="compute",
        account="cfd_project"
    )
    
    # Weak scaling: problem size grows with processor count
    test.set_scaling(
        scaling_type="weak",
        max_nodes=256,
        initial_procs=(4, 4, 4),        # 64 processes initially
        initial_domain=(50.0, 50.0, 50.0),
        initial_cells=(256, 256, 256)    # Grows with domain
    )
    
    return test


def create_gpu_test():
    """
    Create a GPU-accelerated test.
    
    Returns:
        Test: Configured test for GPU execution
    """
    
    test = Test(
        name="cfd_gpu_accelerated",
        input_file=Path("./inputs/cfd_gpu.inp"),
        command=["./cfd_solver_gpu"]
    )
    
    # GPU-specific modules
    test.set_modules([
        "gcc/11.2.0",
        "cuda/11.7",
        "openmpi/4.1.1+cuda"
    ])
    
    # GPU environment
    test.set_env({
        "CUDA_VISIBLE_DEVICES": "0,1,2,3",
        "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
        "NCCL_DEBUG": "INFO"
    })
    
    test.set_backend(scheduler="slurm", launcher="srun")
    
    # GPU resources
    test.set_resources(
        max_nodes=16,
        procs_per_node=4,      # One rank per GPU
        gpus_per_node=4,       # 4 GPUs per node
        time_limit="02:00:00",
        partition="gpu",
        account="gpu_project"
    )
    
    test.set_scaling(
        scaling_type="strong",
        max_nodes=16,
        initial_procs=(2, 2, 1),
        initial_domain=(200.0, 200.0, 100.0),
        initial_cells=(1024, 1024, 512)
    )
    
    return test


# Export the default test
# This is what scaletest.py will find when it loads this file
test = create_cfd_strong_scaling_test()


# Alternative: Export multiple tests
# Users can select which one to run
tests = {
    'strong': create_cfd_strong_scaling_test(),
    'weak': create_cfd_weak_scaling_test(),
    'gpu': create_gpu_test()
}


if __name__ == '__main__':
    """
    When run directly, print test configuration.
    """
    print("=" * 70)
    print("HPC-ScaleTest: Complete Example")
    print("=" * 70)
    print()
    
    # Show strong scaling configuration
    strong_test = create_cfd_strong_scaling_test()
    config = strong_test.to_config()
    
    print("Strong Scaling Test Configuration:")
    print(f"  Name: {config.test_name}")
    print(f"  Input: {config.input_file}")
    print(f"  Command: {' '.join(config.benchmark_command)}")
    print(f"  Scaling: {config.scaling.scaling_type.name}")
    print(f"  Max nodes: {config.scaling.max_nodes}")
    print(f"  Initial procs: {config.scaling.initial_xproc} x "
          f"{config.scaling.initial_yproc} x {config.scaling.initial_zproc}")
    print(f"  Problem size: {config.scaling.initial_nxc}^3 cells")
    print()
    
    print("To run this test:")
    print("  Local:   python scaletest.py run --test complete_example.py --backend local --max-nodes 1")
    print("  Cluster: python scaletest.py run --test complete_example.py --max-nodes 128")
    print()
    print("To validate:")
    print("  python scaletest.py validate --test complete_example.py")
    print()
