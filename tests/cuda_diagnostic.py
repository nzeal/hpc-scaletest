#!/usr/bin/env python3
"""
CUDA Configuration Diagnostic Tool

This script checks your CUDA installation and verifies it's properly
configured for use with the HPC-ScaleTest CMake backend.

Usage:
    python3 cuda_diagnostic.py
"""

import os
import sys
from pathlib import Path
import subprocess

class Colors:
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    END = '\033[0m'
    BOLD = '\033[1m'

def print_header(text):
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{text:^70}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.END}\n")

def print_success(text):
    print(f"{Colors.GREEN}✓ {text}{Colors.END}")

def print_warning(text):
    print(f"{Colors.YELLOW}⚠ {text}{Colors.END}")

def print_error(text):
    print(f"{Colors.RED}✗ {text}{Colors.END}")

def print_info(text):
    print(f"  {text}")

def check_env_variable(var_name):
    """Check if an environment variable is set."""
    value = os.environ.get(var_name)
    if value:
        print_success(f"${var_name} = {value}")
        return value
    else:
        print_warning(f"${var_name} not set")
        return None

def check_path_exists(path, description):
    """Check if a path exists."""
    if path and os.path.exists(path):
        print_success(f"{description}: {path}")
        return True
    else:
        print_error(f"{description} not found: {path}")
        return False

def find_cuda_libraries(cuda_root):
    """Find CUDA libraries in the installation."""
    print("\n📚 Searching for CUDA libraries...")
    
    lib_dirs = []
    for subdir in ['lib64', 'lib', 'lib/x86_64-linux-gnu']:
        lib_path = os.path.join(cuda_root, subdir)
        if os.path.exists(lib_path):
            lib_dirs.append(lib_path)
            print_success(f"Found lib directory: {lib_path}")
    
    if not lib_dirs:
        print_error("No library directories found!")
        return False
    
    # Check for critical libraries
    print("\n🔍 Checking for critical CUDA libraries...")
    critical_libs = {
        'libcudart.so': 'CUDA Runtime (shared)',
        'libcudart_static.a': 'CUDA Runtime (static)',
        'libcudadevrt.a': 'CUDA Device Runtime',
        'libcuda.so': 'CUDA Driver',
    }
    
    found_libs = {}
    for lib_dir in lib_dirs:
        for lib_file, description in critical_libs.items():
            lib_path = os.path.join(lib_dir, lib_file)
            if os.path.exists(lib_path):
                found_libs[lib_file] = lib_path
                print_success(f"{description}: {lib_path}")
    
    # Check if critical libraries are found
    missing = []
    for lib_file, description in critical_libs.items():
        if lib_file not in found_libs:
            print_error(f"{description} not found ({lib_file})")
            missing.append(lib_file)
    
    if 'libcudart_static.a' in missing or 'libcudadevrt.a' in missing:
        print_error("\nCRITICAL: Missing libraries needed for static linking!")
        print_info("This will cause 'cannot find -lcudart_static' or 'cannot find -lcudadevrt' errors")
        return False
    
    return True

def check_nvcc(cuda_root):
    """Check if nvcc is available and working."""
    print("\n🔧 Checking CUDA compiler (nvcc)...")
    
    nvcc_path = os.path.join(cuda_root, 'bin', 'nvcc')
    if not os.path.exists(nvcc_path):
        print_error(f"nvcc not found at: {nvcc_path}")
        return False
    
    print_success(f"nvcc found: {nvcc_path}")
    
    # Try to run nvcc --version
    try:
        result = subprocess.run(
            [nvcc_path, '--version'],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            # Extract version from output
            for line in result.stdout.split('\n'):
                if 'release' in line.lower():
                    print_success(f"nvcc version: {line.strip()}")
                    break
            return True
        else:
            print_error(f"nvcc failed to run (exit code {result.returncode})")
            return False
    except Exception as e:
        print_error(f"Error running nvcc: {e}")
        return False

def check_cmake_available():
    """Check if CMake is available."""
    print("\n🏗️  Checking CMake availability...")
    
    try:
        result = subprocess.run(
            ['cmake', '--version'],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            version_line = result.stdout.split('\n')[0]
            print_success(f"CMake available: {version_line}")
            return True
        else:
            print_error("CMake not working properly")
            return False
    except FileNotFoundError:
        print_error("CMake not found in PATH")
        print_info("Install CMake or load the CMake module")
        return False
    except Exception as e:
        print_error(f"Error checking CMake: {e}")
        return False

def generate_cmake_flags(cuda_root, lib_dirs):
    """Generate CMake flags that would be used."""
    print("\n⚙️  CMake flags that would be generated:")
    print_info("(These flags will be automatically set by the fixed cmake.py)")
    print()
    
    flags = {
        'CMAKE_CUDA_COMPILER': os.path.join(cuda_root, 'bin', 'nvcc'),
        'CMAKE_CUDA_HOST_COMPILER': 'g++',
        'CMAKE_CUDA_COMPILER_FORCED': 'ON',
        'CUDA_TOOLKIT_ROOT_DIR': cuda_root,
        'CMAKE_PREFIX_PATH': cuda_root,
        'CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES': ';'.join(lib_dirs),
    }
    
    # Generate linker flags
    link_flags = []
    for lib_dir in lib_dirs:
        link_flags.append(f"-L{lib_dir}")
    link_flags.extend(["-lcudart", "-lcudadevrt"])
    linker_flags_str = ' '.join(link_flags)
    
    flags['CMAKE_EXE_LINKER_FLAGS'] = linker_flags_str
    flags['CMAKE_SHARED_LINKER_FLAGS'] = linker_flags_str
    
    max_key_len = max(len(k) for k in flags.keys())
    for key, value in flags.items():
        print(f"  {key:{max_key_len}} = {value}")

def find_cuda_root():
    """Find CUDA root directory."""
    print("\n🔍 Searching for CUDA installation...")
    
    # Check environment variables
    for var in ['CUDA_HOME', 'CUDA_ROOT', 'CUDA_PATH', 'CUDADIR']:
        cuda_path = os.environ.get(var)
        if cuda_path and os.path.exists(cuda_path):
            print_success(f"Found CUDA via ${var}: {cuda_path}")
            return cuda_path
    
    # Check common locations
    common_paths = [
        '/usr/local/cuda',
        '/opt/cuda',
        '/usr/cuda',
        '/opt/nvidia/hpc_sdk/Linux_x86_64/cuda',
    ]
    
    for path in common_paths:
        if os.path.exists(path):
            print_success(f"Found CUDA at common location: {path}")
            return path
    
    print_error("CUDA installation not found!")
    print_info("Set CUDA_HOME or CUDA_ROOT environment variable")
    print_info("Or install CUDA in a standard location")
    return None

def main():
    print_header("CUDA Configuration Diagnostic Tool")
    
    print("This tool checks your CUDA installation and verifies it's ready")
    print("for building CUDA applications with the HPC-ScaleTest framework.")
    
    all_checks_passed = True
    
    # Check environment variables
    print_header("Environment Variables")
    cuda_from_env = None
    for var in ['CUDA_HOME', 'CUDA_ROOT', 'CUDA_PATH', 'CUDADIR']:
        value = check_env_variable(var)
        if value and not cuda_from_env:
            cuda_from_env = value
    
    if not cuda_from_env:
        print_warning("No CUDA environment variables set (will use default paths)")
    
    # Find CUDA root
    print_header("CUDA Installation")
    cuda_root = find_cuda_root()
    
    if not cuda_root:
        print_error("\n❌ CRITICAL: CUDA installation not found!")
        print_info("\nACTION REQUIRED:")
        print_info("1. Install CUDA toolkit")
        print_info("2. Set CUDA_HOME environment variable")
        print_info("3. Or load CUDA module: module load cuda")
        return 1
    
    # Check nvcc
    if not check_nvcc(cuda_root):
        print_error("\n❌ CRITICAL: CUDA compiler (nvcc) not working!")
        all_checks_passed = False
    
    # Find and validate libraries
    print_header("CUDA Libraries")
    lib_dirs = []
    for subdir in ['lib64', 'lib', 'lib/x86_64-linux-gnu']:
        lib_path = os.path.join(cuda_root, subdir)
        if os.path.exists(lib_path):
            lib_dirs.append(lib_path)
    
    if not find_cuda_libraries(cuda_root):
        print_error("\n❌ CRITICAL: Required CUDA libraries missing!")
        all_checks_passed = False
    
    # Check CMake
    if not check_cmake_available():
        print_warning("\n⚠️  CMake not available (you may need to load module)")
        print_info("Run: module load cmake")
    
    # Generate example flags
    if lib_dirs:
        generate_cmake_flags(cuda_root, lib_dirs)
    
    # Summary
    print_header("Summary")
    
    if all_checks_passed:
        print_success("✅ All checks passed!")
        print_info("\nYour CUDA installation is properly configured.")
        print_info("The cmake.py backend should work correctly with your setup.")
        print_info("\nNext steps:")
        print_info("1. Ensure the fixed cmake.py is in backends/builds/")
        print_info("2. Run your build with HPC-ScaleTest")
        print_info("3. Monitor build logs for CUDA configuration messages")
        return 0
    else:
        print_error("❌ Some checks failed!")
        print_info("\nYour CUDA installation has issues that need to be resolved.")
        print_info("Review the errors above and take corrective action.")
        print_info("\nCommon solutions:")
        print_info("1. Load CUDA module: module load cuda/12.6")
        print_info("2. Set environment: export CUDA_HOME=/opt/cuda")
        print_info("3. Reinstall CUDA if libraries are missing")
        return 1

if __name__ == '__main__':
    sys.exit(main())
