#!/usr/bin/env python3
"""
Test script to verify the binary detection flow works correctly.

This simulates what happens when:
1. User provides a source directory
2. Framework compiles the code
3. Framework detects the binary name
4. Job script uses the correct binary path
"""

import sys
import os
import tempfile
import shutil
from pathlib import Path

# Add the project to path
sys.path.insert(0, str(Path(__file__).parent))

from backends.builds.cmake import CMakeBackend


def test_binary_detection():
    """Test that CMakeBackend can detect binaries correctly."""
    
    print("=" * 60)
    print("Testing Binary Detection")
    print("=" * 60)
    
    # Create a mock build directory with a fake binary
    with tempfile.TemporaryDirectory() as tmpdir:
        build_dir = Path(tmpdir) / "build"
        build_dir.mkdir()
        
        # Create a fake ELF binary (just the header is enough for detection)
        fake_binary = build_dir / "iPIC3D"
        # ELF magic number: 0x7f 'E' 'L' 'F'
        with open(fake_binary, 'wb') as f:
            f.write(b'\x7fELF' + b'\x00' * 100)
        
        # Make it executable
        fake_binary.chmod(0o755)
        
        # Also create some non-binary files that should be ignored
        (build_dir / "Makefile").write_text("# Makefile\n")
        (build_dir / "CMakeCache.txt").write_text("# CMake cache\n")
        (build_dir / "libtest.so").write_bytes(b'\x7fELF' + b'\x00' * 100)
        (build_dir / "libtest.so").chmod(0o755)
        
        # Create source dir for hints
        source_dir = Path(tmpdir) / "iPIC3D-CPU-NS"
        source_dir.mkdir()
        
        # Test the detection
        backend = CMakeBackend()
        
        print(f"\nBuild directory: {build_dir}")
        print(f"Source directory: {source_dir}")
        print(f"Contents of build dir:")
        for item in build_dir.iterdir():
            print(f"  - {item.name}")
        
        # Test find_executable
        detected = backend.find_executable(build_dir, source_dir)
        
        if detected:
            print(f"\n✓ SUCCESS: Detected binary: {detected.name}")
            print(f"  Full path: {detected}")
            
            # Verify it's the right one (not libtest.so)
            if detected.name == "iPIC3D":
                print(f"  ✓ Correctly identified iPIC3D (not library)")
                return True
            else:
                print(f"  ✗ Wrong binary detected: {detected.name}")
                return False
        else:
            print(f"\n✗ FAILED: No binary detected")
            return False


def test_binary_scoring():
    """Test that binary scoring works for similar names."""
    
    print("\n" + "=" * 60)
    print("Testing Binary Scoring")
    print("=" * 60)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        build_dir = Path(tmpdir)
        
        # Create multiple potential binaries
        binaries = ["iPIC3D", "iPIC3D_test", "benchmark", "main"]
        
        for name in binaries:
            binary = build_dir / name
            with open(binary, 'wb') as f:
                f.write(b'\x7fELF' + b'\x00' * 100)
            binary.chmod(0o755)
        
        source_dir = Path("/fake/path/iPIC3D-CPU-NS")
        
        backend = CMakeBackend()
        detected = backend.find_executable(build_dir, source_dir)
        
        if detected:
            print(f"\n✓ Detected: {detected.name}")
            if detected.name == "iPIC3D":
                print(f"  ✓ Correctly chose iPIC3D over test/benchmark")
                return True
            else:
                print(f"  ⚠ Chose {detected.name} instead of iPIC3D")
                return False
        else:
            print(f"\n✗ No binary detected")
            return False


def test_placeholder_replacement():
    """Test that placeholder gets replaced with actual path."""
    
    print("\n" + "=" * 60)
    print("Testing Placeholder Replacement Logic")
    print("=" * 60)
    
    # Simulate what runner.py does
    command = ["__BINARY_PLACEHOLDER__", "input.txt"]
    binary_path = Path("/workspace/iPIC3D-CPU-NS/build/iPIC3D")
    
    print(f"Original command: {command}")
    print(f"Detected binary: {binary_path}")
    
    # This is what runner.py does
    for i, arg in enumerate(command):
        if arg == "__BINARY_PLACEHOLDER__" or "PLACEHOLDER" in arg:
            command[i] = str(binary_path)
            break
    
    print(f"Updated command: {command}")
    
    if command[0] == str(binary_path):
        print(f"\n✓ SUCCESS: Placeholder correctly replaced")
        return True
    else:
        print(f"\n✗ FAILED: Placeholder not replaced")
        return False


def test_slurm_parsing():
    """Test that slurm.py would correctly parse the binary path."""
    
    print("\n" + "=" * 60)
    print("Testing SLURM Binary Parsing")
    print("=" * 60)
    
    # Simulate what slurm.py does
    command = ["/workspace/iPIC3D-CPU-NS/build/iPIC3D", "os-stdin"]
    
    binary_idx = None
    for idx, arg in enumerate(command):
        if arg.startswith('-') or arg.startswith('$'):
            continue
        if '/' in arg or (arg and arg[0].isalpha()):
            binary_idx = idx
            break
    
    if binary_idx is not None:
        full_binary_path = command[binary_idx]
        binary_path_obj = Path(full_binary_path)
        binary_dir = str(binary_path_obj.parent)
        binary_name = binary_path_obj.name
        
        print(f"Full path: {full_binary_path}")
        print(f"Binary dir: {binary_dir}")
        print(f"Binary name: {binary_name}")
        
        if binary_name == "iPIC3D":
            print(f"\n✓ SUCCESS: Would generate correct job script")
            print(f"  export BINARY={binary_dir}")
            print(f"  srun $BINARY/{binary_name} os-stdin")
            return True
    
    print(f"\n✗ FAILED: Could not parse binary path")
    return False


if __name__ == "__main__":
    results = []
    
    results.append(("Binary Detection", test_binary_detection()))
    results.append(("Binary Scoring", test_binary_scoring()))
    results.append(("Placeholder Replacement", test_placeholder_replacement()))
    results.append(("SLURM Parsing", test_slurm_parsing()))
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    all_passed = True
    for name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status}: {name}")
        if not passed:
            all_passed = False
    
    print()
    if all_passed:
        print("All tests passed!")
        sys.exit(0)
    else:
        print("Some tests failed!")
        sys.exit(1)
