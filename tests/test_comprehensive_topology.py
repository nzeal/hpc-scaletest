#!/usr/bin/env python3
"""
tests/test_comprehensive_topology.py - Comprehensive Test Suite

This test suite verifies:
1. Hardware topology detection
2. MPI configuration derivation
3. SLURM job script generation
4. Correct separation of SLURM allocation vs MPI execution

Key Test Cases:
- Leonardo Booster: 32 cores, 4 GPUs
- LUMI-G: 128 cores, 8 GPUs
- CPU-only partition: 64 cores, 0 GPUs

Author: HPC-ScaleTest Contributors
"""

import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def test_hardware_topology_leonardo():
    """Test hardware topology for Leonardo Booster."""
    print("\n" + "=" * 70)
    print(" TEST: Hardware Topology - Leonardo Booster")
    print("=" * 70)
    
    from core.hardware import HardwareTopology, GPUVendor
    
    topo = HardwareTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        gpu_vendor=GPUVendor.NVIDIA,
        gpu_model="a100",
        detection_method="test"
    )
    
    print(f"\n  Input: 32 CPU cores, 4 GPUs")
    print(f"  Derived: {topo.mpi_ranks_per_node} ranks/node, {topo.cores_per_rank} cores/rank")
    
    assert topo.mpi_ranks_per_node == 4, f"Expected 4, got {topo.mpi_ranks_per_node}"
    assert topo.cores_per_rank == 8, f"Expected 8, got {topo.cores_per_rank}"
    assert topo.is_gpu_partition() == True
    
    print("  ✓ Derived values correct")
    print("  ✓ Leonardo Booster test PASSED")


def test_hardware_topology_lumi():
    """Test hardware topology for LUMI-G."""
    print("\n" + "=" * 70)
    print(" TEST: Hardware Topology - LUMI-G")
    print("=" * 70)
    
    from core.hardware import HardwareTopology, GPUVendor
    
    topo = HardwareTopology(
        partition="standard-g",
        cpu_cores_per_node=128,
        gpus_per_node=8,
        gpu_vendor=GPUVendor.AMD,
        gpu_model="mi250x",
        detection_method="test"
    )
    
    print(f"\n  Input: 128 CPU cores, 8 GPUs")
    print(f"  Derived: {topo.mpi_ranks_per_node} ranks/node, {topo.cores_per_rank} cores/rank")
    
    assert topo.mpi_ranks_per_node == 8, f"Expected 8, got {topo.mpi_ranks_per_node}"
    assert topo.cores_per_rank == 16, f"Expected 16, got {topo.cores_per_rank}"
    
    print("  ✓ LUMI-G test PASSED")


def test_hardware_topology_cpu_only():
    """Test hardware topology for CPU-only partition."""
    print("\n" + "=" * 70)
    print(" TEST: Hardware Topology - CPU Only")
    print("=" * 70)
    
    from core.hardware import HardwareTopology, GPUVendor
    
    topo = HardwareTopology(
        partition="cpu_partition",
        cpu_cores_per_node=64,
        gpus_per_node=0,
        detection_method="test"
    )
    
    print(f"\n  Input: 64 CPU cores, 0 GPUs")
    print(f"  Derived: {topo.mpi_ranks_per_node} ranks/node, {topo.cores_per_rank} cores/rank")
    
    # CPU-only: 1 rank per core
    assert topo.mpi_ranks_per_node == 64
    assert topo.cores_per_rank == 1
    assert topo.is_gpu_partition() == False
    
    print("  ✓ CPU-only test PASSED")


def test_mpi_configuration():
    """Test MPI configuration generation."""
    print("\n" + "=" * 70)
    print(" TEST: MPI Configuration")
    print("=" * 70)
    
    from core.hardware import HardwareTopology, MPIConfiguration, GPUVendor
    
    topo = HardwareTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        detection_method="test"
    )
    
    print("\n  Testing different node counts:")
    
    for N in [1, 2, 4, 8]:
        cfg = MPIConfiguration(num_nodes=N, topology=topo)
        args = cfg.get_mpirun_args()
        cmd = ' '.join(['mpirun'] + args)
        
        expected_np = N * 4
        print(f"    {N} node(s): mpirun -np {cfg.total_mpi_ranks} --map-by ppr:{cfg.mpi_ranks_per_node}:node:PE={cfg.cores_per_rank}")
        
        assert cfg.total_mpi_ranks == expected_np, f"Expected {expected_np}, got {cfg.total_mpi_ranks}"
        assert cfg.mpi_ranks_per_node == 4
        assert cfg.cores_per_rank == 8
        assert f'-np {expected_np}' in cmd
        assert 'ppr:4:node:PE=8' in cmd
    
    print("  ✓ MPI configuration test PASSED")


def test_slurm_script_generation():
    """Test SLURM job script generation with correct resource separation."""
    print("\n" + "=" * 70)
    print(" TEST: SLURM Job Script Generation")
    print("=" * 70)
    
    from core.slurm_script import SlurmJobConfig, generate_slurm_job_script
    
    # Leonardo Booster configuration
    config = SlurmJobConfig(
        job_name="test_job",
        num_nodes=4,
        partition="boost_usr_prod",
        cpus_per_node=32,       # Full node allocation
        gpus_per_node=4,
        mpi_ranks_per_node=4,   # Actual MPI ranks
        cores_per_rank=8,
        account="my_account",
        time_limit="02:00:00",
    )
    
    script = generate_slurm_job_script(
        config,
        "$BINARY/iPIC3D",
        ["os-stdin"],
        ["module load cuda/12.0"]
    )
    
    print("\n  Verifying SLURM resource allocation:")
    
    # CRITICAL: Check that --ntasks-per-node uses FULL CPU count
    assert "#SBATCH --ntasks-per-node=32" in script, \
        "FAILED: --ntasks-per-node should be 32 (full node), not 4 (MPI ranks)"
    print("  ✓ --ntasks-per-node=32 (full node allocation)")
    
    assert "#SBATCH --gres=gpu:4" in script
    print("  ✓ --gres=gpu:4 (full GPU allocation)")
    
    print("\n  Verifying MPI execution parameters:")
    
    # Check MPI command uses actual ranks
    assert "mpirun -np 16" in script, "MPI should use 16 ranks (4 nodes × 4 ranks/node)"
    print("  ✓ mpirun -np 16 (actual MPI ranks)")
    
    assert "ppr:4:node:PE=8" in script, "MPI mapping should be 4 ranks/node with 8 cores each"
    print("  ✓ --map-by ppr:4:node:PE=8")
    
    assert "./bind.sh" in script, "GPU binding script should be included"
    print("  ✓ GPU binding script included")
    
    assert "CUDA_VISIBLE_DEVICES" in script, "Script should set CUDA_VISIBLE_DEVICES"
    print("  ✓ CUDA_VISIBLE_DEVICES binding")
    
    print("\n  ✓ SLURM job script test PASSED")


def test_expected_commands_leonardo():
    """
    Test that generated commands match the expected format for Leonardo Booster.
    
    Expected behavior (from requirements):
    - 1 node: mpirun -np 4 --map-by ppr:4:node:PE=8
    - 2 nodes: mpirun -np 8 --map-by ppr:4:node:PE=8
    - 4 nodes: mpirun -np 16 --map-by ppr:4:node:PE=8
    """
    print("\n" + "=" * 70)
    print(" TEST: Expected MPI Commands (Leonardo)")
    print("=" * 70)
    
    from core.hardware import HardwareTopology, MPIConfiguration
    
    topo = HardwareTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        detection_method="test"
    )
    
    expected = {
        1: "mpirun -np 4 --map-by ppr:4:node:PE=8 --report-bindings",
        2: "mpirun -np 8 --map-by ppr:4:node:PE=8 --report-bindings",
        4: "mpirun -np 16 --map-by ppr:4:node:PE=8 --report-bindings",
    }
    
    print("\n  Comparing generated vs expected commands:")
    
    for N, expected_cmd in expected.items():
        cfg = MPIConfiguration(num_nodes=N, topology=topo)
        actual_cmd = 'mpirun ' + cfg.get_mpirun_command_string()
        
        print(f"\n    {N} node(s):")
        print(f"      Expected: {expected_cmd}")
        print(f"      Actual:   {actual_cmd}")
        
        # Check key components
        assert f'-np {N * 4}' in actual_cmd, f"Wrong -np for {N} nodes"
        assert 'ppr:4:node:PE=8' in actual_cmd, f"Wrong ppr mapping for {N} nodes"
        
    print("\n  ✓ Expected commands test PASSED")


def test_portability_8gpu_system():
    """
    Test portability to systems with 8 GPUs per node (like DGX or LUMI).
    
    This must work WITHOUT code changes.
    """
    print("\n" + "=" * 70)
    print(" TEST: Portability - 8 GPU System")
    print("=" * 70)
    
    from core.hardware import HardwareTopology, MPIConfiguration
    from core.slurm_script import SlurmJobConfig, generate_slurm_job_script
    
    # DGX-style or LUMI configuration
    topo = HardwareTopology(
        partition="gpu",
        cpu_cores_per_node=128,
        gpus_per_node=8,
        detection_method="test"
    )
    
    print(f"\n  System: 128 cores, 8 GPUs per node")
    print(f"  Derived: {topo.mpi_ranks_per_node} ranks/node, {topo.cores_per_rank} cores/rank")
    
    assert topo.mpi_ranks_per_node == 8
    assert topo.cores_per_rank == 16
    
    # Test job config
    config = SlurmJobConfig(
        job_name="test_8gpu",
        num_nodes=4,
        partition="gpu",
        cpus_per_node=128,
        gpus_per_node=8,
        mpi_ranks_per_node=8,
        cores_per_rank=16,
    )
    
    script = generate_slurm_job_script(config, "./app", ["input.dat"])
    
    # Verify
    assert "#SBATCH --ntasks-per-node=128" in script
    print("  ✓ --ntasks-per-node=128")
    
    assert "#SBATCH --gres=gpu:8" in script
    print("  ✓ --gres=gpu:8")
    
    assert "mpirun -np 32" in script  # 4 × 8
    print("  ✓ mpirun -np 32")
    
    assert "ppr:8:node:PE=16" in script
    print("  ✓ --map-by ppr:8:node:PE=16")
    
    print("\n  ✓ 8-GPU portability test PASSED")


def main():
    """Run all tests."""
    print("=" * 70)
    print(" HPC-ScaleTest - Comprehensive Topology & Job Script Tests")
    print("=" * 70)
    
    try:
        test_hardware_topology_leonardo()
        test_hardware_topology_lumi()
        test_hardware_topology_cpu_only()
        test_mpi_configuration()
        test_slurm_script_generation()
        test_expected_commands_leonardo()
        test_portability_8gpu_system()
        
        print("\n" + "=" * 70)
        print(" ALL TESTS PASSED")
        print("=" * 70)
        
        print("\n  Summary of verified behaviors:")
        print("    ✓ Hardware topology detection")
        print("    ✓ MPI rank/core derivation (ranks = GPUs, cores = CPUs/GPUs)")
        print("    ✓ SLURM: --ntasks-per-node uses FULL CPU count")
        print("    ✓ SLURM: --gres=gpu:N uses full GPU allocation")
        print("    ✓ MPI: -np = nodes × ranks_per_node")
        print("    ✓ MPI: --map-by ppr:R:node:PE=C syntax")
        print("    ✓ GPU binding via CUDA_VISIBLE_DEVICES")
        print("    ✓ Portable across different GPU counts (4, 8 GPUs)")
        
        return 0
        
    except AssertionError as e:
        print(f"\n  FAILED: {e}")
        return 1
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
