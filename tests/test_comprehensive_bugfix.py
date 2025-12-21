#!/usr/bin/env python3
"""
tests/test_comprehensive_bugfix.py

Comprehensive test suite verifying all bug fixes and requirements:

1. Bug Identification and Fixes
   - Scaling behavior bugs
   - Job script generation bugs
   - MPI launch configuration bugs
   - GPU binding bugs

2. Correct Slurm and MPI Integration
   - Job script correctness
   - bind.sh generation
   - GPU binding via CUDA_VISIBLE_DEVICES

3. No Hardcoded Hardware Assumptions
   - CPUs per node detected at runtime
   - GPUs per node detected at runtime
   - Works across heterogeneous systems

4. Correct MPI Mapping Logic
   - Ranks per node = GPUs (for GPU jobs)
   - Cores per rank = CPUs / GPUs
   - Consistent between SLURM and mpirun

5. Node Sequence Generation
   - Powers of 2 up to max_nodes
   - No hardcoded [1, 2, 4]

6. Output and Data Management
   - run_data/ directory
   - Per-run JSON files
   - CSV output
   - Aggregated reports

Author: HPC-ScaleTest Contributors
"""

import sys
import os
import json
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


def print_section(title):
    print(f"\n{'=' * 70}")
    print(f" {title}")
    print('=' * 70)


def print_test(name):
    print(f"\n  [{name}]")


def print_pass(msg):
    print(f"    ✓ {msg}")


def print_fail(msg):
    print(f"    ✗ {msg}")
    raise AssertionError(msg)


# =============================================================================
# 1. BUG IDENTIFICATION AND FIXES
# =============================================================================

def test_scaling_behavior_bugs():
    """Verify scaling behavior bugs are fixed."""
    print_section("1. Bug Identification and Fixes")
    print_test("Scaling Behavior")
    
    from core.unified_execution import (
        HardwareTopology, MPIConfiguration, generate_node_sequence
    )
    
    # Bug: Hardcoded node sequences
    # Fix: Dynamic generation based on max_nodes
    seq = generate_node_sequence(max_nodes=8)
    assert seq == [1, 2, 4, 8], f"Node sequence bug: {seq}"
    print_pass("Node sequence dynamically generated (not hardcoded)")
    
    # Bug: Incorrect MPI rank calculation
    # Fix: ranks_per_node = gpus_per_node (for GPU jobs)
    topo = HardwareTopology(
        partition="test",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        detection_method="test"
    )
    assert topo.mpi_ranks_per_node == 4, "MPI ranks should equal GPUs"
    print_pass("MPI ranks = GPUs per node (for GPU jobs)")
    
    # Bug: Total ranks don't scale correctly
    # Fix: total_ranks = nodes × ranks_per_node
    for nodes in [1, 2, 4, 8]:
        config = MPIConfiguration(num_nodes=nodes, topology=topo)
        expected = nodes * 4
        assert config.total_mpi_ranks == expected, f"Total ranks bug: {config.total_mpi_ranks} != {expected}"
    print_pass("Total MPI ranks scale correctly with node count")


def test_job_script_generation_bugs():
    """Verify job script generation bugs are fixed."""
    print_test("Job Script Generation")
    
    from core.unified_execution import (
        HardwareTopology, MPIConfiguration, SlurmConfiguration,
        generate_job_script
    )
    
    topo = HardwareTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        detection_method="test"
    )
    
    mpi = MPIConfiguration(num_nodes=4, topology=topo)
    slurm = SlurmConfiguration(
        job_name="test",
        num_nodes=4,
        partition="boost_usr_prod",
        topology=topo
    )
    
    script = generate_job_script(
        slurm_config=slurm,
        mpi_config=mpi,
        binary_dir="/path/to/build",
        executable="iPIC3D",
        args=["os-stdin"]
    )
    
    # Bug: BINARY=${PWD}
    # Fix: BINARY must be absolute path to build directory
    assert "export BINARY=/path/to/build" in script, "BINARY bug: should be absolute path"
    assert "export BINARY=${PWD}" not in script, "BINARY bug: should not be ${PWD}"
    print_pass("BINARY points to absolute build directory (not ${PWD})")
    
    # Bug: Executable is bind.sh
    # Fix: Executable must be the actual binary
    assert "$BINARY/iPIC3D" in script, "Executable bug: should be actual binary"
    print_pass("Executable is actual binary ($BINARY/iPIC3D)")
    
    # Bug: mpirun missing mapping
    # Fix: mpirun must include ppr:R:node:PE=C
    assert "ppr:4:node:PE=8" in script, "mpirun mapping bug"
    print_pass("mpirun includes ppr:R:node:PE=C mapping")


def test_mpi_launch_configuration_bugs():
    """Verify MPI launch configuration bugs are fixed."""
    print_test("MPI Launch Configuration")
    
    from core.unified_execution import HardwareTopology, MPIConfiguration
    
    # Test Leonardo Booster configuration
    topo = HardwareTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        detection_method="test"
    )
    
    config = MPIConfiguration(num_nodes=4, topology=topo)
    cmd = config.get_mpirun_command_string()
    
    # Bug: Wrong -np value
    # Fix: -np = nodes × ranks_per_node
    assert "-np 16" in cmd, f"-np bug: {cmd}"
    print_pass("-np = 16 (4 nodes × 4 ranks/node)")
    
    # Bug: Wrong mapping syntax
    # Fix: --map-by ppr:R:node:PE=C
    assert "--map-by ppr:4:node:PE=8" in cmd, f"mapping bug: {cmd}"
    print_pass("--map-by ppr:4:node:PE=8 (correct syntax)")
    
    # Bug: Missing --report-bindings
    # Fix: Always include for debugging
    assert "--report-bindings" in cmd, f"report-bindings missing: {cmd}"
    print_pass("--report-bindings included")


def test_gpu_binding_bugs():
    """Verify GPU binding bugs are fixed."""
    print_test("GPU Binding")
    
    from core.unified_execution import (
        GPU_BIND_SCRIPT, generate_job_script,
        HardwareTopology, MPIConfiguration, SlurmConfiguration
    )
    
    # Bug: bind.sh not generated
    # Fix: bind.sh generated inline with correct content
    assert "OMPI_COMM_WORLD_LOCAL_RANK" in GPU_BIND_SCRIPT, "bind.sh missing OMPI rank"
    assert "CUDA_VISIBLE_DEVICES=$LOCAL_RANK" in GPU_BIND_SCRIPT, "bind.sh missing CUDA binding"
    print_pass("bind.sh uses OMPI_COMM_WORLD_LOCAL_RANK for GPU binding")
    
    # Bug: bind.sh has hardcoded values
    # Fix: Uses environment variables with fallbacks
    assert "MV2_COMM_WORLD_LOCAL_RANK" in GPU_BIND_SCRIPT, "bind.sh missing MVAPICH2 fallback"
    assert "SLURM_LOCALID" in GPU_BIND_SCRIPT, "bind.sh missing SLURM fallback"
    print_pass("bind.sh has MPI implementation fallbacks")
    
    # Bug: bind.sh not included in job script
    # Fix: Generated inline and wraps executable
    topo = HardwareTopology(
        partition="test",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        detection_method="test"
    )
    mpi = MPIConfiguration(num_nodes=4, topology=topo)
    slurm = SlurmConfiguration(
        job_name="test", num_nodes=4,
        partition="test", topology=topo
    )
    script = generate_job_script(slurm, mpi, "/build", "app", [])
    
    assert "cat > bind.sh << 'BIND_EOF'" in script, "bind.sh not generated in script"
    assert "./bind.sh $BINARY/app" in script, "bind.sh not wrapping executable"
    print_pass("bind.sh generated inline and wraps executable")


# =============================================================================
# 2. CORRECT SLURM AND MPI INTEGRATION
# =============================================================================

def test_slurm_mpi_integration():
    """Verify SLURM and MPI integration is correct."""
    print_section("2. Correct Slurm and MPI Integration")
    print_test("SLURM Directives")
    
    from core.unified_execution import (
        HardwareTopology, SlurmConfiguration
    )
    
    topo = HardwareTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        detection_method="test"
    )
    
    slurm = SlurmConfiguration(
        job_name="test",
        num_nodes=4,
        partition="boost_usr_prod",
        topology=topo
    )
    
    directives = slurm.get_directives()
    directives_str = '\n'.join(directives)
    
    # Bug: --ntasks-per-node uses MPI ranks
    # Fix: --ntasks-per-node uses FULL CPU count
    assert "--ntasks-per-node=32" in directives_str, "SLURM bug: should use full CPU count"
    print_pass("--ntasks-per-node=32 (full CPU count, not MPI ranks)")
    
    # Bug: Missing GPU directive
    # Fix: --gres=gpu:N included for GPU partitions
    assert "--gres=gpu:4" in directives_str, "SLURM bug: missing GPU directive"
    print_pass("--gres=gpu:4 (GPU allocation)")
    
    # Verify CPU-only partition doesn't request GPUs
    cpu_topo = HardwareTopology(
        partition="cpu_partition",
        cpu_cores_per_node=64,
        gpus_per_node=0,
        detection_method="test"
    )
    cpu_slurm = SlurmConfiguration(
        job_name="cpu_test",
        num_nodes=4,
        partition="cpu_partition",
        topology=cpu_topo
    )
    cpu_directives = '\n'.join(cpu_slurm.get_directives())
    
    assert "--gres=gpu" not in cpu_directives, "CPU-only should not request GPUs"
    assert "--ntasks-per-node=64" in cpu_directives
    print_pass("CPU-only partition: no GPU directive, full CPU allocation")


# =============================================================================
# 3. NO HARDCODED HARDWARE ASSUMPTIONS
# =============================================================================

def test_no_hardcoded_values():
    """Verify no hardcoded hardware assumptions."""
    print_section("3. No Hardcoded Hardware Assumptions")
    print_test("Dynamic Hardware Detection")
    
    from core.unified_execution import HardwareTopology, MPIConfiguration
    
    # Test various configurations - all should work without code changes
    test_systems = [
        # (name, cpus, gpus, expected_ranks, expected_cores)
        ("Leonardo Booster", 32, 4, 4, 8),
        ("LUMI-G", 128, 8, 8, 16),
        ("DGX A100", 128, 8, 8, 16),
        ("Custom 6-GPU", 48, 6, 6, 8),
        ("Custom 16-GPU", 128, 16, 16, 8),
        ("CPU-only 64 core", 64, 0, 64, 1),
        ("CPU-only 128 core", 128, 0, 128, 1),
    ]
    
    for name, cpus, gpus, exp_ranks, exp_cores in test_systems:
        topo = HardwareTopology(
            partition="test",
            cpu_cores_per_node=cpus,
            gpus_per_node=gpus,
            detection_method="test"
        )
        
        assert topo.mpi_ranks_per_node == exp_ranks, f"{name}: wrong ranks"
        assert topo.cores_per_rank == exp_cores, f"{name}: wrong cores/rank"
        print_pass(f"{name}: {cpus} CPUs, {gpus} GPUs → {exp_ranks} ranks, {exp_cores} cores/rank")
    
    print_test("MPI Scaling")
    
    # Verify MPI scaling works for all configurations
    for name, cpus, gpus, exp_ranks, exp_cores in test_systems[:3]:  # First 3
        topo = HardwareTopology(
            partition="test",
            cpu_cores_per_node=cpus,
            gpus_per_node=gpus,
            detection_method="test"
        )
        
        for nodes in [1, 2, 4, 8]:
            config = MPIConfiguration(num_nodes=nodes, topology=topo)
            expected_total = nodes * exp_ranks
            assert config.total_mpi_ranks == expected_total
        
        print_pass(f"{name}: MPI ranks scale correctly with nodes")


# =============================================================================
# 4. CORRECT MPI MAPPING LOGIC
# =============================================================================

def test_mpi_mapping_logic():
    """Verify MPI mapping logic is correct."""
    print_section("4. Correct MPI Mapping Logic")
    print_test("Leonardo Booster Example")
    
    from core.unified_execution import HardwareTopology, MPIConfiguration
    
    # Leonardo Booster: 32 CPUs, 4 GPUs
    topo = HardwareTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        detection_method="test"
    )
    
    # Verify derived configuration
    assert topo.mpi_ranks_per_node == 4, "Ranks should equal GPUs"
    assert topo.cores_per_rank == 8, "Cores/rank = CPUs/GPUs = 32/4 = 8"
    print_pass("Derived: 4 ranks/node, 8 cores/rank")
    
    # Verify expected MPI commands
    expected_commands = {
        1: "mpirun -np 4 --map-by ppr:4:node:PE=8 --report-bindings",
        2: "mpirun -np 8 --map-by ppr:4:node:PE=8 --report-bindings",
        4: "mpirun -np 16 --map-by ppr:4:node:PE=8 --report-bindings",
    }
    
    for nodes, expected in expected_commands.items():
        config = MPIConfiguration(num_nodes=nodes, topology=topo)
        actual = config.get_mpirun_command_string()
        assert actual == expected, f"Command mismatch for {nodes} nodes"
        print_pass(f"{nodes} node(s): {actual}")
    
    print_test("SLURM/MPI Consistency")
    
    from core.unified_execution import SlurmConfiguration
    
    # Verify SLURM and MPI are consistent
    slurm = SlurmConfiguration(
        job_name="test",
        num_nodes=4,
        partition="boost_usr_prod",
        topology=topo
    )
    
    directives = '\n'.join(slurm.get_directives())
    
    # SLURM requests full node resources
    assert "--ntasks-per-node=32" in directives
    print_pass("SLURM: --ntasks-per-node=32 (full node)")
    
    # MPI uses derived configuration
    config = MPIConfiguration(num_nodes=4, topology=topo)
    assert config.total_mpi_ranks == 16
    assert config.mpi_ranks_per_node == 4
    print_pass("MPI: 16 total ranks (4 per node)")
    
    # GPU binding matches
    assert config.total_gpus == 16
    print_pass("GPUs: 16 total (4 per node × 4 nodes)")


# =============================================================================
# 5. NODE SEQUENCE GENERATION
# =============================================================================

def test_node_sequence_generation():
    """Verify node sequence generation."""
    print_section("5. Node Sequence Generation")
    print_test("Powers of 2")
    
    from core.unified_execution import generate_node_sequence
    
    # Test various max_nodes values
    test_cases = [
        (1, [1]),
        (2, [1, 2]),
        (4, [1, 2, 4]),
        (8, [1, 2, 4, 8]),
        (16, [1, 2, 4, 8, 16]),
        (32, [1, 2, 4, 8, 16, 32]),
    ]
    
    for max_nodes, expected in test_cases:
        actual = generate_node_sequence(max_nodes)
        assert actual == expected, f"max_nodes={max_nodes}: {actual} != {expected}"
        print_pass(f"max_nodes={max_nodes}: {actual}")
    
    print_test("Non-Power-of-2 Max")
    
    # Non-power-of-2 includes the max value
    seq = generate_node_sequence(max_nodes=6)
    assert seq == [1, 2, 4, 6], f"Expected [1, 2, 4, 6], got {seq}"
    print_pass("max_nodes=6: [1, 2, 4, 6] (includes 6)")
    
    seq = generate_node_sequence(max_nodes=12)
    assert seq == [1, 2, 4, 8, 12], f"Expected [1, 2, 4, 8, 12], got {seq}"
    print_pass("max_nodes=12: [1, 2, 4, 8, 12] (includes 12)")


# =============================================================================
# 6. OUTPUT AND DATA MANAGEMENT
# =============================================================================

def test_output_data_management():
    """Verify output and data management."""
    print_section("6. Output and Data Management")
    print_test("Directory Structure")
    
    from core.unified_execution import (
        IncrementalResultWriter, HardwareTopology, MPIConfiguration
    )
    
    topo = HardwareTopology(
        partition="boost_usr_prod",
        cpu_cores_per_node=32,
        gpus_per_node=4,
        detection_method="test"
    )
    
    with tempfile.TemporaryDirectory() as tmpdir:
        writer = IncrementalResultWriter(
            output_dir=Path(tmpdir),
            scaling_type="weak",
            partition="boost_usr_prod",
            test_name="test"
        )
        
        # Verify run_data directory created
        data_dir = Path(tmpdir) / "run_data"
        assert data_dir.exists(), "run_data/ not created"
        print_pass("run_data/ directory created")
        
        print_test("Per-Run JSON Files")
        
        # Simulate runs
        for nodes in [1, 2, 4]:
            config = MPIConfiguration(num_nodes=nodes, topology=topo)
            result = writer.create_pending_result(
                node_count=nodes,
                mpi_config=config,
                topology=topo,
                job_directory=f"{tmpdir}/node{nodes}"
            )
            
            # Verify file written immediately
            expected_file = data_dir / f"weak_scaling_report_node{nodes}.json"
            assert expected_file.exists(), f"Per-run JSON not written: {expected_file}"
            
            # Verify content
            with open(expected_file) as f:
                data = json.load(f)
            assert data['node_count'] == nodes
            assert data['gpus_per_node'] == 4
            assert data['total_gpus'] == nodes * 4
            assert 'timestamp' in data
        
        print_pass("Per-run JSON files written immediately")
        
        print_test("Run Data Content")
        
        # Update with completion data
        for nodes in [1, 2, 4]:
            writer.update_run_status(
                node_count=nodes,
                status="completed",
                wall_time=100.0 + nodes * 5,
                exit_code=0
            )
        
        # Verify content updated
        with open(data_dir / "weak_scaling_report_node1.json") as f:
            data = json.load(f)
        
        assert data['status'] == 'completed'
        assert data['wall_time_seconds'] > 0
        assert data['total_mpi_ranks'] == 4
        assert data['mpi_ranks_per_node'] == 4
        assert data['cores_per_rank'] == 8
        print_pass("Run data includes all required fields")
        
        print_test("Summary and Reports")
        
        # Generate summary and reports
        writer.write_summary()
        writer.write_csv_report()
        writer.write_aggregated_report()
        
        # Verify all files exist
        assert (data_dir / "summary.json").exists()
        assert (data_dir / "scaling_report.csv").exists()
        assert (data_dir / "scaling_report.json").exists()
        print_pass("summary.json created")
        print_pass("scaling_report.csv created")
        print_pass("scaling_report.json created")
        
        # Verify summary content
        with open(data_dir / "summary.json") as f:
            summary = json.load(f)
        
        assert summary['completed_runs'] == 3
        assert summary['total_runs'] == 3
        assert summary['node_sequence'] == [1, 2, 4]
        print_pass("Summary contains correct run counts")
        
        # Verify CSV content
        with open(data_dir / "scaling_report.csv") as f:
            csv_lines = f.read().strip().split('\n')
        
        assert len(csv_lines) == 4  # Header + 3 data rows
        assert 'nodes,gpus_per_node,total_gpus' in csv_lines[0]
        print_pass("CSV contains header and data rows")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    print("\n" + "=" * 70)
    print(" HPC-ScaleTest - Comprehensive Bug Fix Verification")
    print("=" * 70)
    
    tests = [
        test_scaling_behavior_bugs,
        test_job_script_generation_bugs,
        test_mpi_launch_configuration_bugs,
        test_gpu_binding_bugs,
        test_slurm_mpi_integration,
        test_no_hardcoded_values,
        test_mpi_mapping_logic,
        test_node_sequence_generation,
        test_output_data_management,
    ]
    
    passed = 0
    failed = 0
    
    for test_func in tests:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"\n  ✗ TEST FAILED: {e}")
            failed += 1
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 70)
    print(f" RESULTS: {passed}/{len(tests)} test suites passed")
    print("=" * 70)
    
    if failed == 0:
        print("\n  ✓ ALL BUG FIXES VERIFIED")
        print("\n  Summary of Fixed Bugs:")
        print("    1. Scaling behavior bugs (node sequence, rank calculation)")
        print("    2. Job script generation (BINARY path, executable, mapping)")
        print("    3. MPI launch configuration (-np, --map-by, bindings)")
        print("    4. GPU binding (bind.sh, CUDA_VISIBLE_DEVICES)")
        print("    5. SLURM/MPI integration (ntasks-per-node, consistency)")
        print("    6. Hardcoded values eliminated (all values derived)")
        print("    7. Node sequence generation (dynamic, powers of 2)")
        print("    8. Output management (incremental, per-run JSON, CSV)")
        print()
    else:
        sys.exit(1)
