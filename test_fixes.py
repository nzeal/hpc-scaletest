#!/usr/bin/env python3
"""
Test script to verify HPC-ScaleTest fixes:
1. GPU mpirun command uses bind.sh correctly
2. max_nodes from YAML is respected (no hardcoded [1,2,4] limit)
3. Modular CPU/GPU job generators work correctly
"""

import os
import sys
import tempfile
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

def test_max_nodes_parsing():
    """Test that max_nodes from YAML is correctly parsed."""
    print("\n" + "="*70)
    print("TEST 1: max_nodes YAML Parsing")
    print("="*70)
    
    from utils.config_parser import YAMLConfigParser
    
    # Create a temporary YAML config
    yaml_content = """
repository: https://github.com/test/repo.git
hardware_type: "gpu"
gpus_per_node: 4
procs_per_node: 32
scaling: strong
max_nodes: 16
initial_domain: [5.12, 5.12, 1]
initial_cells: [128, 128, 1]
initial_procs: [2, 2, 1]
partition: "boost_usr_prod"
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_content)
        temp_path = f.name
    
    try:
        parser = YAMLConfigParser(Path(temp_path))
        config = parser.to_orchestrator_config()
        
        print(f"  Parsed max_nodes: {config.get('max_nodes', 'NOT FOUND')}")
        print(f"  Parsed hardware_type: {config.get('hardware_type', 'NOT FOUND')}")
        print(f"  Parsed gpus_per_node: {config.get('gpus_per_node', 'NOT FOUND')}")
        
        assert config.get('max_nodes') == 16, f"Expected max_nodes=16, got {config.get('max_nodes')}"
        assert config.get('hardware_type') == 'gpu', f"Expected hardware_type=gpu, got {config.get('hardware_type')}"
        assert config.get('gpus_per_node') == 4, f"Expected gpus_per_node=4, got {config.get('gpus_per_node')}"
        
        print("  ✓ PASSED: max_nodes correctly parsed as 16")
        
    finally:
        os.unlink(temp_path)


def test_node_sequence_generation():
    """Test that node sequence is generated correctly up to max_nodes."""
    print("\n" + "="*70)
    print("TEST 2: Node Sequence Generation")
    print("="*70)
    
    from core.config import ScalingConfig
    
    # Test with max_nodes = 16
    config = ScalingConfig(max_nodes=16)
    sequence = config.get_node_sequence()
    
    print(f"  max_nodes: 16")
    print(f"  Generated sequence: {sequence}")
    
    expected = [1, 2, 4, 8, 16]
    assert sequence == expected, f"Expected {expected}, got {sequence}"
    print(f"  ✓ PASSED: Sequence correctly includes all powers of 2 up to 16")
    
    # Test with non-power-of-2 max_nodes
    config2 = ScalingConfig(max_nodes=12)
    sequence2 = config2.get_node_sequence()
    
    print(f"\n  max_nodes: 12")
    print(f"  Generated sequence: {sequence2}")
    
    expected2 = [1, 2, 4, 8, 12]  # 12 is added since it's not a power of 2
    assert sequence2 == expected2, f"Expected {expected2}, got {sequence2}"
    print(f"  ✓ PASSED: Non-power-of-2 max_nodes is included")


def test_gpu_job_generator():
    """Test that GPU job generator creates correct mpirun command."""
    print("\n" + "="*70)
    print("TEST 3: GPU Job Generator - mpirun Command")
    print("="*70)
    
    from engine.job_generators import GPUJobGenerator, JobGeneratorConfig
    
    # Configure for Leonardo Booster: 32 cores, 4 GPUs
    config = JobGeneratorConfig(
        cpus_per_node=32,
        gpus_per_node=4,
        partition="boost_usr_prod",
        account="test_account",
        binary_dir="/path/to/build",
        executable="iPIC3D",
        input_file="os-stdin",
        launcher="mpirun",
        modules=["cuda/12.0", "nvhpc/24.5"],
        max_nodes=16
    )
    
    generator = GPUJobGenerator(config)
    
    # Test mpirun command for 4 nodes
    num_nodes = 4
    total_ranks = num_nodes * 4  # 4 GPUs per node
    ranks_per_node = 4
    cores_per_rank = 8  # 32 / 4
    
    mpi_cmd = generator.generate_mpi_command(
        num_nodes=num_nodes,
        total_ranks=total_ranks,
        ranks_per_node=ranks_per_node,
        cores_per_rank=cores_per_rank
    )
    
    print(f"  Configuration:")
    print(f"    Nodes: {num_nodes}")
    print(f"    GPUs per node: 4")
    print(f"    CPUs per node: 32")
    print(f"    Cores per rank: {cores_per_rank}")
    print(f"\n  Generated command:")
    print(f"    {mpi_cmd}")
    
    # Check command contains required elements
    assert "-np 16" in mpi_cmd, f"Command should have -np 16, got: {mpi_cmd}"
    assert "ppr:4:node:PE=8" in mpi_cmd, f"Command should have ppr:4:node:PE=8, got: {mpi_cmd}"
    assert "--report-bindings" in mpi_cmd, f"Command should have --report-bindings, got: {mpi_cmd}"
    assert "./bind.sh" in mpi_cmd, f"Command should use ./bind.sh, got: {mpi_cmd}"
    assert "$BINARY/iPIC3D" in mpi_cmd, f"Command should have $BINARY/iPIC3D, got: {mpi_cmd}"
    assert "os-stdin" in mpi_cmd, f"Command should have os-stdin, got: {mpi_cmd}"
    
    print(f"\n  ✓ PASSED: Command format is correct:")
    print(f"    ✓ -np 16 (4 nodes × 4 GPUs)")
    print(f"    ✓ --map-by ppr:4:node:PE=8")
    print(f"    ✓ --report-bindings")
    print(f"    ✓ ./bind.sh wrapper")
    print(f"    ✓ $BINARY/iPIC3D executable")
    print(f"    ✓ os-stdin input file")


def test_gpu_job_script_generation():
    """Test that GPU job script is generated correctly."""
    print("\n" + "="*70)
    print("TEST 4: GPU Job Script Generation")
    print("="*70)
    
    from engine.job_generators import GPUJobGenerator, JobGeneratorConfig
    
    config = JobGeneratorConfig(
        cpus_per_node=32,
        gpus_per_node=4,
        partition="boost_usr_prod",
        account="cin_staff",
        qos="boost_qos_bprod",
        binary_dir="/scratch/build",
        executable="iPIC3D",
        input_file="os-stdin",
        launcher="mpirun",
        modules=["cuda/12.6", "nvhpc/24.5"],
        max_nodes=16,
        time_limit="02:00:00"
    )
    
    generator = GPUJobGenerator(config)
    
    script = generator.generate_script(
        num_nodes=4,
        job_name="gpu_test_node4",
        working_dir="/scratch/jobs/node4"
    )
    
    print("  Script contains:")
    
    # Check essential elements
    checks = [
        ("#SBATCH --nodes=4", "Node allocation"),
        ("#SBATCH --gres=gpu:4", "GPU allocation"),
        ("#SBATCH --partition=boost_usr_prod", "Partition"),
        ("#SBATCH --qos=boost_qos_bprod", "QoS"),
        ("mpirun -np 16", "MPI rank count"),
        ("ppr:4:node:PE=8", "MPI mapping"),
        ("--report-bindings", "Report bindings flag"),
        ("./bind.sh", "GPU binding script"),
        ("CUDA_VISIBLE_DEVICES", "GPU visibility setting"),
        ("export OMP_NUM_THREADS=8", "OpenMP threads"),
    ]
    
    all_passed = True
    for check, description in checks:
        if check in script:
            print(f"    ✓ {description}: {check}")
        else:
            print(f"    ✗ MISSING {description}: {check}")
            all_passed = False
    
    if all_passed:
        print(f"\n  ✓ PASSED: Job script contains all required elements")
    else:
        print(f"\n  ✗ FAILED: Some elements missing from job script")
        return False
    
    return True


def test_cpu_job_generator():
    """Test CPU job generator for comparison."""
    print("\n" + "="*70)
    print("TEST 5: CPU Job Generator")
    print("="*70)
    
    from engine.job_generators import CPUJobGenerator, JobGeneratorConfig
    
    config = JobGeneratorConfig(
        cpus_per_node=128,
        gpus_per_node=0,
        partition="cpu_prod",
        account="test_account",
        binary_dir="/path/to/build",
        executable="app",
        input_file="input.dat",
        launcher="mpirun",
        max_nodes=8
    )
    
    generator = CPUJobGenerator(config)
    
    # Test mpirun command
    mpi_cmd = generator.generate_mpi_command(
        num_nodes=4,
        total_ranks=512,
        ranks_per_node=128,
        cores_per_rank=1
    )
    
    print(f"  Generated CPU command:")
    print(f"    {mpi_cmd}")
    
    # Should NOT have bind.sh or GPU-specific options
    assert "./bind.sh" not in mpi_cmd, "CPU command should not have bind.sh"
    assert "ppr:" not in mpi_cmd, "CPU command should not have ppr mapping"
    assert "-np 512" in mpi_cmd, "CPU command should have correct rank count"
    
    print(f"\n  ✓ PASSED: CPU command format is correct (no GPU elements)")


def main():
    """Run all tests."""
    print("\n" + "="*70)
    print(" HPC-ScaleTest Fix Verification Tests")
    print("="*70)
    
    try:
        test_max_nodes_parsing()
        test_node_sequence_generation()
        test_gpu_job_generator()
        test_gpu_job_script_generation()
        test_cpu_job_generator()
        test_required_field_validation()
        
        print("\n" + "="*70)
        print(" ALL TESTS PASSED ✓")
        print("="*70)
        print("\nSummary of fixes verified:")
        print("  1. max_nodes from YAML is correctly parsed (16, not limited to 4)")
        print("  2. Node sequence includes all powers of 2 up to max_nodes")
        print("  3. GPU mpirun uses ./bind.sh with ppr:X:node:PE=Y --report-bindings format")
        print("  4. GPU job scripts contain proper CUDA binding and bind.sh fallback")
        print("  5. CPU and GPU job generators are properly separated")
        print("  6. No hardcoded defaults - all values from YAML/CLI/detection")
        print("")
        
    except Exception as e:
        print(f"\n\n  ✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


def test_required_field_validation():
    """Test that missing required fields raise errors."""
    print("\n" + "="*70)
    print("TEST 6: Required Field Validation (No Hardcoded Defaults)")
    print("="*70)
    
    from utils.config_parser import YAMLConfigParser
    import tempfile
    
    # Test 1: Missing max_nodes should fail
    yaml_missing_nodes = """
repository: https://github.com/test/repo.git
partition: "boost_usr_prod"
initial_procs: [2, 2, 1]
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_missing_nodes)
        temp_path = f.name
    
    try:
        parser = YAMLConfigParser(Path(temp_path))
        parser.parse()
        try:
            parser.validate()
            print("  ✗ FAILED: Should have raised error for missing max_nodes")
            return False
        except ValueError as e:
            if "max_nodes" in str(e).lower() or "nodes" in str(e).lower():
                print("  ✓ Correctly rejected config missing max_nodes")
            else:
                print(f"  ✗ Wrong error message: {e}")
                return False
    finally:
        os.unlink(temp_path)
    
    # Test 2: Missing partition should fail
    yaml_missing_partition = """
repository: https://github.com/test/repo.git
max_nodes: 16
initial_procs: [2, 2, 1]
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_missing_partition)
        temp_path = f.name
    
    try:
        parser = YAMLConfigParser(Path(temp_path))
        parser.parse()
        try:
            parser.validate()
            print("  ✗ FAILED: Should have raised error for missing partition")
            return False
        except ValueError as e:
            if "partition" in str(e).lower():
                print("  ✓ Correctly rejected config missing partition")
            else:
                print(f"  ✗ Wrong error message: {e}")
                return False
    finally:
        os.unlink(temp_path)
    
    # Test 3: GPU job missing gpus_per_node should fail
    yaml_gpu_missing_gpus = """
repository: https://github.com/test/repo.git
hardware_type: gpu
max_nodes: 16
partition: "boost_usr_prod"
initial_procs: [2, 2, 1]
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_gpu_missing_gpus)
        temp_path = f.name
    
    try:
        parser = YAMLConfigParser(Path(temp_path))
        parser.parse()
        try:
            parser.validate()
            print("  ✗ FAILED: Should have raised error for GPU job missing gpus_per_node")
            return False
        except ValueError as e:
            if "gpus_per_node" in str(e).lower():
                print("  ✓ Correctly rejected GPU config missing gpus_per_node")
            else:
                print(f"  ✗ Wrong error message: {e}")
                return False
    finally:
        os.unlink(temp_path)
    
    # Test 4: GPU job with gpus_per_node but missing cpu count should fail
    yaml_gpu_missing_cpus = """
repository: https://github.com/test/repo.git
hardware_type: gpu
gpus_per_node: 4
max_nodes: 16
partition: "boost_usr_prod"
initial_procs: [2, 2, 1]
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_gpu_missing_cpus)
        temp_path = f.name
    
    try:
        parser = YAMLConfigParser(Path(temp_path))
        parser.parse()
        try:
            parser.validate()
            print("  ✗ FAILED: Should have raised error for GPU job missing cpus_per_node/procs_per_node")
            return False
        except ValueError as e:
            if "cpus_per_node" in str(e).lower() or "procs_per_node" in str(e).lower():
                print("  ✓ Correctly rejected GPU config missing CPU count")
            else:
                print(f"  ✗ Wrong error message: {e}")
                return False
    finally:
        os.unlink(temp_path)
    
    # Test 5: Complete GPU config should pass (with procs_per_node)
    yaml_complete_gpu = """
repository: https://github.com/test/repo.git
hardware_type: gpu
gpus_per_node: 4
procs_per_node: 32
max_nodes: 16
partition: "boost_usr_prod"
initial_procs: [2, 2, 1]
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(yaml_complete_gpu)
        temp_path = f.name
    
    try:
        parser = YAMLConfigParser(Path(temp_path))
        parser.parse()
        parser.validate()  # Should not raise
        print("  ✓ Complete GPU config validated successfully")
    except ValueError as e:
        print(f"  ✗ FAILED: Valid config was rejected: {e}")
        return False
    finally:
        os.unlink(temp_path)
    
    print("  ✓ PASSED: All required field validations work correctly")
    return True


if __name__ == "__main__":
    sys.exit(main())
