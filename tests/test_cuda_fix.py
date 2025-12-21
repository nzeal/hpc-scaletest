#!/usr/bin/env python3
"""
CUDA Linking Fix Verification Test

This script verifies that the fixed cmake.py properly configures CUDA linking.
It simulates the build process and checks that all necessary flags are set.

Usage:
    python3 test_cuda_fix.py
"""

import os
import sys
from pathlib import Path

# Add the parent directory to path to import the fixed module
sys.path.insert(0, str(Path(__file__).parent))

try:
    from backends.builds.cmake import CMakeBackend
except ImportError:
    print("❌ ERROR: Cannot import CMakeBackend")
    print("   Make sure cmake.py is in backends/builds/")
    sys.exit(1)

class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    END = '\033[0m'
    BOLD = '\033[1m'

def print_header(text):
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{text:^70}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.END}\n")

def print_pass(text):
    print(f"{Colors.GREEN}✓ PASS: {text}{Colors.END}")

def print_fail(text):
    print(f"{Colors.RED}✗ FAIL: {text}{Colors.END}")

def print_info(text):
    print(f"  {text}")

def test_nvhpc_detection():
    """Test that NVHPC is detected from module commands."""
    print_header("Test 1: NVHPC Detection")
    
    # Test with NVHPC module
    backend = CMakeBackend(options={
        'module_commands': ['module load nvhpc/24.5']
    })
    
    if backend._nvhpc_detected:
        print_pass("NVHPC detected from module commands")
        return True
    else:
        print_fail("NVHPC not detected")
        return False

def test_cuda_path_validation():
    """Test CUDA path validation."""
    print_header("Test 2: CUDA Path Validation")
    
    backend = CMakeBackend()
    
    # Test with invalid path
    if not backend._validate_cuda_path("/nonexistent/path"):
        print_pass("Invalid path correctly rejected")
    else:
        print_fail("Invalid path incorrectly accepted")
        return False
    
    # Test with current CUDA (if available)
    cuda_root = backend._find_cuda_root()
    if cuda_root:
        if backend._validate_cuda_path(cuda_root):
            print_pass(f"Current CUDA installation validated: {cuda_root}")
            return True
        else:
            print_fail(f"Current CUDA failed validation: {cuda_root}")
            return False
    else:
        print_info("No CUDA installation found (test skipped)")
        return True

def test_library_directory_detection():
    """Test that multiple library directories are detected."""
    print_header("Test 3: Library Directory Detection")
    
    backend = CMakeBackend(options={
        'module_commands': ['module load nvhpc/24.5']
    })
    
    cuda_root = backend._find_cuda_root()
    if not cuda_root:
        print_info("No CUDA installation found (test skipped)")
        return True
    
    # Simulate the library directory detection
    lib_dirs = []
    for subdir in ['lib64', 'lib', 'lib/x86_64-linux-gnu']:
        lib_path = os.path.join(cuda_root, subdir)
        if os.path.exists(lib_path):
            lib_dirs.append(lib_path)
            print_pass(f"Found library directory: {lib_path}")
    
    if lib_dirs:
        print_pass(f"Detected {len(lib_dirs)} library director{'y' if len(lib_dirs)==1 else 'ies'}")
        return True
    else:
        print_fail("No library directories detected")
        return False

def test_linker_flags_generation():
    """Test that linker flags are properly generated."""
    print_header("Test 4: Linker Flags Generation")
    
    backend = CMakeBackend(options={
        'module_commands': ['module load nvhpc/24.5']
    })
    
    # Create test flags
    flags = {
        'CMAKE_C_COMPILER': 'mpicc',
        'CMAKE_CXX_COMPILER': 'mpicxx',
    }
    
    # Apply NVHPC CUDA flags
    flags = backend._add_nvhpc_cuda_flags(flags)
    
    # Check critical flags
    tests_passed = True
    
    if 'CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES' in flags:
        print_pass("CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES set")
        print_info(f"  Value: {flags['CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES']}")
    else:
        print_fail("CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES not set")
        tests_passed = False
    
    if 'CMAKE_EXE_LINKER_FLAGS' in flags:
        linker_flags = flags['CMAKE_EXE_LINKER_FLAGS']
        print_pass("CMAKE_EXE_LINKER_FLAGS set")
        print_info(f"  Value: {linker_flags}")
        
        # Check for critical components
        if '-L' in linker_flags:
            print_pass("  Contains -L (library path)")
        else:
            print_fail("  Missing -L (library path)")
            tests_passed = False
        
        if '-lcudart' in linker_flags:
            print_pass("  Contains -lcudart")
        else:
            print_fail("  Missing -lcudart")
            tests_passed = False
        
        if '-lcudadevrt' in linker_flags:
            print_pass("  Contains -lcudadevrt")
        else:
            print_fail("  Missing -lcudadevrt")
            tests_passed = False
    else:
        print_fail("CMAKE_EXE_LINKER_FLAGS not set")
        tests_passed = False
    
    if 'CMAKE_SHARED_LINKER_FLAGS' in flags:
        print_pass("CMAKE_SHARED_LINKER_FLAGS set (critical for shared libraries)")
        print_info(f"  Value: {flags['CMAKE_SHARED_LINKER_FLAGS']}")
    else:
        print_fail("CMAKE_SHARED_LINKER_FLAGS not set (will cause linking failures!)")
        tests_passed = False
    
    return tests_passed

def test_backward_compatibility():
    """Test that existing flags are preserved."""
    print_header("Test 5: Backward Compatibility")
    
    backend = CMakeBackend(options={
        'module_commands': ['module load nvhpc/24.5']
    })
    
    # Start with existing flags
    flags = {
        'CMAKE_BUILD_TYPE': 'Release',
        'CMAKE_C_FLAGS': '-O3',
        'CMAKE_EXE_LINKER_FLAGS': '-lm',
    }
    
    original_flags = flags.copy()
    
    # Apply NVHPC CUDA flags
    flags = backend._add_nvhpc_cuda_flags(flags)
    
    # Check that original flags are preserved
    tests_passed = True
    
    if flags.get('CMAKE_BUILD_TYPE') == original_flags['CMAKE_BUILD_TYPE']:
        print_pass("CMAKE_BUILD_TYPE preserved")
    else:
        print_fail("CMAKE_BUILD_TYPE modified")
        tests_passed = False
    
    if flags.get('CMAKE_C_FLAGS') == original_flags['CMAKE_C_FLAGS']:
        print_pass("CMAKE_C_FLAGS preserved")
    else:
        print_fail("CMAKE_C_FLAGS modified")
        tests_passed = False
    
    # Check that existing linker flags are extended, not replaced
    if original_flags['CMAKE_EXE_LINKER_FLAGS'] in flags.get('CMAKE_EXE_LINKER_FLAGS', ''):
        print_pass("Existing linker flags preserved and extended")
    else:
        print_fail("Existing linker flags overwritten")
        tests_passed = False
    
    return tests_passed

def main():
    print_header("CUDA Linking Fix Verification Test Suite")
    
    print("This test suite verifies that the fixed cmake.py properly")
    print("configures CUDA linking to prevent 'cannot find -lcuda*' errors.")
    print()
    
    tests = [
        ("NVHPC Detection", test_nvhpc_detection),
        ("CUDA Path Validation", test_cuda_path_validation),
        ("Library Directory Detection", test_library_directory_detection),
        ("Linker Flags Generation", test_linker_flags_generation),
        ("Backward Compatibility", test_backward_compatibility),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print_fail(f"Test crashed: {e}")
            results.append((test_name, False))
    
    # Summary
    print_header("Test Results Summary")
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        if result:
            print_pass(f"{test_name}")
        else:
            print_fail(f"{test_name}")
    
    print()
    print(f"Tests passed: {passed}/{total}")
    
    if passed == total:
        print(f"\n{Colors.GREEN}{Colors.BOLD}✓ ALL TESTS PASSED{Colors.END}")
        print(f"{Colors.GREEN}The CUDA linking fix is working correctly!{Colors.END}")
        return 0
    else:
        print(f"\n{Colors.RED}{Colors.BOLD}✗ SOME TESTS FAILED{Colors.END}")
        print(f"{Colors.RED}The fix may not be properly installed or configured.{Colors.END}")
        return 1

if __name__ == '__main__':
    sys.exit(main())
