#!/usr/bin/env python3
"""
System checker for HPC Auto.
Helps diagnose common issues before running tests.
"""

import subprocess
import sys
import os
from pathlib import Path


def check_command(cmd, name):
    """Check if a command is available."""
    try:
        result = subprocess.run(
            [cmd, "--version"],
            capture_output=True,
            text=True,
            timeout=5
        )
        print(f"✓ {name}: Available")
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print(f"✗ {name}: Not found")
        return False


def check_module_system():
    """Check if module system is available."""
    try:
        result = subprocess.run(
            ["module", "avail"],
            capture_output=True,
            text=True,
            timeout=5,
            shell=True
        )
        print(f"✓ Module System: Available")
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print(f"✗ Module System: Not found")
        return False


def check_scheduler():
    """Check which scheduler is available."""
    schedulers = {
        "sbatch": "Slurm",
        "qsub": "PBS/Torque",
    }
    
    found = []
    for cmd, name in schedulers.items():
        try:
            result = subprocess.run(
                [cmd, "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            print(f"✓ Scheduler: {name} (use --scheduler {name.lower()})")
            found.append(name.lower())
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
    
    if not found:
        print(f"✗ Scheduler: None found (use --scheduler local for testing)")
        return False
    
    return True


def check_modules_available(modules):
    """Check if specific modules are available."""
    print(f"\nChecking for common HPC modules:")
    try:
        result = subprocess.run(
            ["module", "avail"],
            capture_output=True,
            text=True,
            timeout=10,
            shell=True
        )
        output = result.stdout + result.stderr
        
        for module in modules:
            base_name = module.split('/')[0]
            if base_name.lower() in output.lower():
                print(f"✓ {base_name}: Available")
            else:
                print(f"✗ {base_name}: Not found")
    except Exception as e:
        print(f"  Could not check modules: {e}")


def suggest_fixes():
    """Suggest fixes for common issues."""
    print("\n" + "="*60)
    print("COMMON SOLUTIONS")
    print("="*60)
    
    print("\n1. If modules are missing:")
    print("   - Use --modules to specify available modules:")
    print("     python hpc_auto.py /path --modules 'gcc/11.2.0,openmpi/4.1.1'")
    print("   - Or skip auto-detection with --no-build")
    
    print("\n2. If build fails:")
    print("   - Check the README has build instructions")
    print("   - Manually specify build system:")
    print("     python hpc_auto.py /path --build-system cmake")
    print("   - Use pre-built executable with --no-build")
    
    print("\n3. If no scheduler is found:")
    print("   - Use local scheduler for testing:")
    print("     python hpc_auto.py /path --scheduler local --nodes 1")
    
    print("\n4. For module issues:")
    print("   - Check available modules: module avail")
    print("   - Load modules before running:")
    print("     module load gcc openmpi")
    print("     python hpc_auto.py /path --nodes 4")


def main():
    print("="*60)
    print("HPC AUTO - SYSTEM CHECK")
    print("="*60)
    print("\nChecking system requirements...\n")
    
    # Check basic commands
    print("Essential Tools:")
    check_command("python", "Python")
    check_command("git", "Git")
    check_command("cmake", "CMake")
    check_command("make", "Make")
    check_command("gcc", "GCC")
    check_command("mpicc", "MPI")
    
    print("\nHPC Environment:")
    check_module_system()
    check_scheduler()
    
    # Check common modules
    common_modules = ["gcc", "openmpi", "mpich", "intel", "cuda", "hdf5"]
    check_modules_available(common_modules)
    
    # Suggest fixes
    suggest_fixes()
    
    print("\n" + "="*60)
    print("System check complete!")
    print("="*60)
    print("\nTo run a test:")
    print("  python hpc_auto.py /path/to/code --nodes 4 --verbose")


if __name__ == '__main__':
    main()
