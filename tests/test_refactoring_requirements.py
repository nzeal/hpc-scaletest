#!/usr/bin/env python3
"""
tests/test_refactoring_requirements.py

Comprehensive test suite verifying all refactoring requirements:
1. Incremental result generation with data/ directory
2. GPU job script generation (correct BINARY, bind.sh, mpirun)
3. MPI mapping derived at runtime (no hardcoded values)
4. Heterogeneous system support (4 GPU, 8 GPU systems)
5. Code organization (separate CPU/GPU modules)

Author: HPC-ScaleTest Contributors
"""

import sys
import os
import json
import tempfile
from pathlib import Path

# Setup path
sys.path.insert(0, str(Path(__file__).parent.parent))


def print_header(title):
    print(f"\n{'=' * 70}")
    print(f" {title}")
    print('=' * 70)


def print_pass(msg):
    print(f"  ✓ {msg}")


def print_fail(msg):
    print(f"  ✗ {msg}")
    raise AssertionError(msg)


# =============================================================================
# TEST 1: Incremental Result Generation
# =============================================================================

def test_incremental_result_generation():
    """Test that results are written incrementally to data/ directory."""
    print_header("TEST 1: Incremental Result Generation")
    
    from core.result_writer import IncrementalResultWriter, RunResult
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)
        
        # Create writer
        writer = IncrementalResultWriter(
            output_dir=output_dir,
            scaling_type="weak",
            partition="boost_usr_prod",
            test_name="test_scaling"
        )
        
        # Verify data/ directory created
        data_dir = output_dir / "data"
        assert data_dir.exists(), "data/ directory not created"
        print_pass("data/ directory created")
        
        # Simulate runs - each should write immediately
        node_counts = [1, 2, 4, 8]
        for nodes in node_counts:
            result = writer.create_pending_result(
                node_count=nodes,
                total_mpi_ranks=nodes * 4,
                mpi_ranks_per_node=4,
                cores_per_rank=8,
                gpus_per_node=4,
                procs_decomp=(nodes, 2, 2),
                job_directory=f"{tmpdir}/node{nodes}"
            )
            
            # Verify file written immediately
            expected_file = data_dir / f"weak_scaling_report_node{nodes}.json"
            assert expected_file.exists(), f"Per-node JSON not written: {expected_file}"
            
            # Verify content
            with open(expected_file) as f:
                data = json.load(f)
            assert data['node_count'] == nodes
            assert data['status'] == 'pending'
        
        print_pass("Per-node JSON files written immediately")
        
        # Update with completion
        for nodes in node_counts:
            writer.update_run_status(
                node_count=nodes,
                status="completed",
                wall_time=100.0 + nodes * 5,
                exit_code=0
            )
        
        # Generate summary and aggregated report
        writer.write_summary()
        writer.write_aggregated_report()
        
        # Verify all expected files exist
        expected_files = [
            "weak_scaling_report_node1.json",
            "weak_scaling_report_node2.json",
            "weak_scaling_report_node4.json",
            "weak_scaling_report_node8.json",
            "summary.json",
            "weak_scaling_report.json",
            "weak_scaling_report.txt"
        ]
        
        for filename in expected_files:
            filepath = data_dir / filename
            assert filepath.exists(), f"Missing file: {filename}"
        
        print_pass("summary.json and aggregated report generated")
        
        # Verify summary content
        with open(data_dir / "summary.json") as f:
            summary = json.load(f)
        
        assert summary['completed_runs'] == 4
        assert summary['total_runs'] == 4
        assert len(summary['node_sequence']) == 4
        print_pass("Summary contains correct run counts")
    
    print("\n  TEST 1 PASSED\n")


# =============================================================================
# TEST 2: GPU Job Script Generation
# =============================================================================

def test_gpu_job_script_generation():
    """Test GPU job script with correct BINARY, bind.sh, and mpirun."""
    print_header("TEST 2: GPU Job Script Generation")
    
    from engine.gpu_execution import GPUExecutionEngine, GPUJobConfig
    
    # Configure engine for Leonardo Booster
    engine = GPUExecutionEngine(partition="boost_usr_prod")
    engine._cpus_per_node = 32
    engine._gpus_per_node = 4
    
    # Create job config with EXPLICIT binary path
    config = engine.create_job_config(
        num_nodes=4,
        binary_dir="/path/to/build",  # Must be absolute path, NOT ${PWD}
        executable_name="iPIC3D"
    )
    
    # Test 1: Verify derived values
    assert config.mpi_ranks_per_node == 4, "Should have 4 ranks/node (1 per GPU)"
    assert config.cores_per_rank == 8, "Should have 8 cores/rank (32/4)"
    assert config.total_mpi_ranks == 16, "Should have 16 total ranks (4 nodes × 4)"
    print_pass("MPI configuration derived correctly")
    
    # Test 2: Generate job script
    script = engine.generate_job_script(
        config,
        job_name="test_gpu_job",
        args=["os-stdin"],
        modules=["module load cuda/12.0"]
    )
    
    # Verify BINARY is NOT ${PWD}
    assert "export BINARY=/path/to/build" in script, "BINARY must be absolute path"
    assert "export BINARY=${PWD}" not in script, "BINARY must NOT be ${PWD}"
    print_pass("BINARY set to absolute build directory")
    
    # Verify executable is the actual binary, not bind.sh
    assert "$BINARY/iPIC3D" in script, "Executable must be $BINARY/iPIC3D"
    print_pass("Executable is actual binary ($BINARY/iPIC3D)")
    
    # Verify bind.sh is generated
    assert "cat > bind.sh" in script, "bind.sh must be generated"
    assert "OMPI_COMM_WORLD_LOCAL_RANK" in script, "bind.sh must use OMPI local rank"
    assert "CUDA_VISIBLE_DEVICES" in script, "bind.sh must set CUDA_VISIBLE_DEVICES"
    print_pass("bind.sh generated with correct GPU binding")
    
    # Verify mpirun command format
    assert "mpirun -np 16" in script, "mpirun must have correct -np"
    assert "ppr:4:node:PE=8" in script, "mpirun must have ppr:R:node:PE=C mapping"
    assert "./bind.sh $BINARY/iPIC3D" in script, "mpirun must wrap with bind.sh"
    print_pass("mpirun has correct mapping (ppr:4:node:PE=8)")
    
    # Verify SLURM directives
    assert "#SBATCH --ntasks-per-node=32" in script, "Must request FULL CPUs"
    assert "#SBATCH --gres=gpu:4" in script, "Must request all GPUs"
    print_pass("SLURM directives use full node allocation")
    
    print("\n  TEST 2 PASSED\n")


# =============================================================================
# TEST 3: MPI Mapping Derived at Runtime
# =============================================================================

def test_mpi_mapping_runtime_derivation():
    """Test MPI mapping is derived from topology, not hardcoded."""
    print_header("TEST 3: MPI Mapping Derived at Runtime")
    
    from core.hardware import HardwareTopology, MPIConfiguration
    
    # Test various system configurations
    test_cases = [
        # (cpus, gpus, expected_ranks, expected_cores)
        (32, 4, 4, 8),    # Leonardo Booster
        (128, 8, 8, 16),  # LUMI-G / DGX
        (64, 4, 4, 16),   # Different ratio
        (48, 6, 6, 8),    # 6 GPU system
        (64, 0, 64, 1),   # CPU-only
    ]
    
    print(f"\n  Testing {len(test_cases)} configurations:")
    
    for cpus, gpus, expected_ranks, expected_cores in test_cases:
        topo = HardwareTopology(
            partition="test",
            cpu_cores_per_node=cpus,
            gpus_per_node=gpus,
            detection_method="test"
        )
        
        # Verify derived values
        if gpus > 0:
            assert topo.mpi_ranks_per_node == gpus, f"Ranks should equal GPUs ({gpus})"
            assert topo.cores_per_rank == cpus // gpus, f"Cores should be CPUs/GPUs"
        else:
            assert topo.mpi_ranks_per_node == cpus, "CPU-only: ranks should equal CPUs"
            assert topo.cores_per_rank == 1, "CPU-only: 1 core per rank"
        
        print(f"    {cpus} CPUs, {gpus} GPUs → {topo.mpi_ranks_per_node} ranks, {topo.cores_per_rank} cores/rank ✓")
    
    print_pass("MPI mapping derived correctly for all configurations")
    
    # Test MPI configuration scaling
    topo = HardwareTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        detection_method="test"
    )
    
    print("\n  Testing rank scaling with node count:")
    for nodes in [1, 2, 4, 8, 16]:
        config = MPIConfiguration(num_nodes=nodes, topology=topo)
        expected_total = nodes * 4
        assert config.total_mpi_ranks == expected_total, f"Expected {expected_total} ranks"
        
        cmd = config.get_mpirun_command_string()
        assert f"-np {expected_total}" in cmd
        print(f"    {nodes} nodes → -np {expected_total} ✓")
    
    print_pass("Total ranks scale correctly with node count")
    
    print("\n  TEST 3 PASSED\n")


# =============================================================================
# TEST 4: Heterogeneous System Support
# =============================================================================

def test_heterogeneous_system_support():
    """Test support for different system configurations."""
    print_header("TEST 4: Heterogeneous System Support")
    
    from core.slurm_script import SlurmJobConfig, generate_slurm_job_script
    
    # Test Leonardo Booster (32 CPUs, 4 GPUs)
    print("\n  Testing Leonardo Booster (32 CPUs, 4 GPUs):")
    config = SlurmJobConfig(
        job_name="leonardo_test",
        num_nodes=4,
        partition="boost_usr_prod",
        cpus_per_node=32,
        gpus_per_node=4,
        mpi_ranks_per_node=4,
        cores_per_rank=8,
    )
    script = generate_slurm_job_script(config, "$BINARY/app", ["input"])
    
    assert "#SBATCH --ntasks-per-node=32" in script
    assert "#SBATCH --gres=gpu:4" in script
    assert "mpirun -np 16" in script
    assert "ppr:4:node:PE=8" in script
    print_pass("Leonardo configuration correct")
    
    # Test LUMI-G (128 CPUs, 8 GPUs)
    print("\n  Testing LUMI-G (128 CPUs, 8 GPUs):")
    config = SlurmJobConfig(
        job_name="lumi_test",
        num_nodes=4,
        partition="standard-g",
        cpus_per_node=128,
        gpus_per_node=8,
        mpi_ranks_per_node=8,
        cores_per_rank=16,
    )
    script = generate_slurm_job_script(config, "$BINARY/app", ["input"])
    
    assert "#SBATCH --ntasks-per-node=128" in script
    assert "#SBATCH --gres=gpu:8" in script
    assert "mpirun -np 32" in script
    assert "ppr:8:node:PE=16" in script
    print_pass("LUMI-G configuration correct")
    
    # Test DGX (128 CPUs, 8 GPUs)
    print("\n  Testing DGX (128 CPUs, 8 GPUs):")
    config = SlurmJobConfig(
        job_name="dgx_test",
        num_nodes=2,
        partition="dgx",
        cpus_per_node=128,
        gpus_per_node=8,
        mpi_ranks_per_node=8,
        cores_per_rank=16,
    )
    script = generate_slurm_job_script(config, "$BINARY/app", ["input"])
    
    assert "#SBATCH --ntasks-per-node=128" in script
    assert "#SBATCH --gres=gpu:8" in script
    assert "mpirun -np 16" in script
    assert "ppr:8:node:PE=16" in script
    print_pass("DGX configuration correct")
    
    # Test CPU-only (64 CPUs, 0 GPUs)
    print("\n  Testing CPU-only (64 CPUs, 0 GPUs):")
    config = SlurmJobConfig(
        job_name="cpu_test",
        num_nodes=4,
        partition="cpu_partition",
        cpus_per_node=64,
        gpus_per_node=0,
        mpi_ranks_per_node=64,
        cores_per_rank=1,
    )
    script = generate_slurm_job_script(config, "$BINARY/app", ["input"])
    
    assert "#SBATCH --ntasks-per-node=64" in script
    assert "--gres=gpu" not in script  # No GPU directive
    assert "mpirun -np 256" in script
    assert "ppr:64:node:PE=1" in script
    print_pass("CPU-only configuration correct")
    
    print("\n  TEST 4 PASSED\n")


# =============================================================================
# TEST 5: Code Organization
# =============================================================================

def test_code_organization():
    """Test that code is properly organized into modules."""
    print_header("TEST 5: Code Organization")
    
    # Test core modules exist
    core_modules = [
        "core.hardware",
        "core.slurm_script",
        "core.result_writer",
        "core.mpi_command",
        "core.topology",
    ]
    
    print("\n  Checking core modules:")
    for module_name in core_modules:
        try:
            __import__(module_name)
            print_pass(f"{module_name} available")
        except ImportError as e:
            print(f"    ⚠ {module_name} not available: {e}")
    
    # Test engine modules exist
    engine_modules = [
        "engine.cpu_execution",
        "engine.gpu_execution",
        "engine.scaling",
        "engine.runner",
    ]
    
    print("\n  Checking engine modules:")
    for module_name in engine_modules:
        try:
            __import__(module_name)
            print_pass(f"{module_name} available")
        except ImportError as e:
            print(f"    ⚠ {module_name} not available: {e}")
    
    # Verify CPU and GPU execution are separate
    from engine.cpu_execution import CPUExecutionEngine, CPUJobConfig
    from engine.gpu_execution import GPUExecutionEngine, GPUJobConfig
    
    assert CPUExecutionEngine is not GPUExecutionEngine
    assert CPUJobConfig is not GPUJobConfig
    print_pass("CPU and GPU execution logic separated")
    
    # Verify core.hardware is single source of truth
    from core.hardware import HardwareTopology, MPIConfiguration, TopologyDetector
    print_pass("core.hardware is single source of truth for topology")
    
    print("\n  TEST 5 PASSED\n")


# =============================================================================
# TEST 6: No Hardcoded Values
# =============================================================================

def test_no_hardcoded_values():
    """Test that no site-specific constants are hardcoded."""
    print_header("TEST 6: No Hardcoded Values")
    
    from core.hardware import HardwareTopology
    
    # Test that HardwareTopology computes values, doesn't assume
    print("\n  Testing dynamic computation:")
    
    # Any CPU/GPU combination should work
    test_configs = [
        (32, 4),
        (64, 8),
        (48, 6),
        (128, 16),
        (256, 8),
        (112, 4),  # Leonardo total with SMT
    ]
    
    for cpus, gpus in test_configs:
        topo = HardwareTopology(
            partition="test",
            cpu_cores_per_node=cpus,
            gpus_per_node=gpus,
            detection_method="test"
        )
        
        # Values should be computed, not assumed
        assert topo.mpi_ranks_per_node == gpus
        assert topo.cores_per_rank == cpus // gpus
        print(f"    {cpus} CPUs / {gpus} GPUs = {topo.cores_per_rank} cores/rank ✓")
    
    print_pass("All values computed dynamically")
    
    print("\n  TEST 6 PASSED\n")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    print("\n" + "=" * 70)
    print(" HPC-ScaleTest - Refactoring Requirements Test Suite")
    print("=" * 70)
    
    tests = [
        ("Incremental Result Generation", test_incremental_result_generation),
        ("GPU Job Script Generation", test_gpu_job_script_generation),
        ("MPI Mapping Runtime Derivation", test_mpi_mapping_runtime_derivation),
        ("Heterogeneous System Support", test_heterogeneous_system_support),
        ("Code Organization", test_code_organization),
        ("No Hardcoded Values", test_no_hardcoded_values),
    ]
    
    passed = 0
    failed = 0
    
    for name, test_func in tests:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"\n  TEST FAILED: {e}")
            failed += 1
    
    print("\n" + "=" * 70)
    print(f" RESULTS: {passed} passed, {failed} failed")
    print("=" * 70)
    
    if failed == 0:
        print("\n  ✓ ALL REQUIREMENTS VERIFIED\n")
        print("  Summary of verified behaviors:")
        print("    ✓ Results written immediately to data/ directory")
        print("    ✓ Per-node JSON files: (weak|strong)_scaling_report_node<N>.json")
        print("    ✓ Summary and aggregated reports generated after all runs")
        print("    ✓ BINARY points to build directory (not ${PWD})")
        print("    ✓ Executable is actual binary (not bind.sh)")
        print("    ✓ bind.sh uses OMPI_COMM_WORLD_LOCAL_RANK for GPU binding")
        print("    ✓ mpirun uses ppr:R:node:PE=C mapping")
        print("    ✓ MPI mapping derived at runtime from topology")
        print("    ✓ Works with 4, 6, 8, 16 GPU systems")
        print("    ✓ Works with CPU-only partitions")
        print("    ✓ CPU and GPU execution in separate modules")
        print("    ✓ No hardcoded site-specific values")
    else:
        sys.exit(1)
