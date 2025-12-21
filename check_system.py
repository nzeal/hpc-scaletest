#!/usr/bin/env python3
"""
System Check Utility for HPC-ScaleTest
Validates system requirements and installation status.
"""

import sys
import subprocess
from pathlib import Path
from typing import Tuple, List


def check_python_version() -> Tuple[bool, str]:
    """Check if Python version meets requirements."""
    version = sys.version_info
    required = (3, 8)
    
    if version >= required:
        return True, f"✓ Python {version.major}.{version.minor}.{version.micro}"
    else:
        return False, f"✗ Python {version.major}.{version.minor} (requires >= 3.8)"


def check_module(module_name: str, required: bool = True) -> Tuple[bool, str]:
    """Check if a Python module is available."""
    try:
        __import__(module_name)
        return True, f"✓ {module_name}"
    except ImportError:
        status = "✗" if required else "○"
        note = " (required)" if required else " (optional)"
        return not required, f"{status} {module_name}{note}"


def check_command(command: str, required: bool = False) -> Tuple[bool, str]:
    """Check if a system command is available."""
    try:
        result = subprocess.run(
            ['which', command], 
            capture_output=True, 
            text=True
        )
        if result.returncode == 0:
            path = result.stdout.strip()
            return True, f"✓ {command} ({path})"
        else:
            status = "✗" if required else "○"
            note = " (required)" if required else " (optional)"
            return not required, f"{status} {command}{note}"
    except Exception:
        return False, f"✗ {command}"


def check_scheduler() -> Tuple[bool, str, str]:
    """Detect available job scheduler."""
    schedulers = {
        'sbatch': 'SLURM',
        'qsub': 'PBS/Torque',
    }
    
    for cmd, name in schedulers.items():
        try:
            result = subprocess.run(
                ['which', cmd],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                return True, name, f"✓ Scheduler: {name}"
        except Exception:
            pass
    
    return True, 'Local', "○ No scheduler detected (local execution only)"


def check_mpi() -> Tuple[bool, str]:
    """Check for MPI installation."""
    mpi_commands = ['mpirun', 'mpiexec', 'srun']
    
    found = []
    for cmd in mpi_commands:
        try:
            result = subprocess.run(
                ['which', cmd],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                found.append(cmd)
        except Exception:
            pass
    
    if found:
        return True, f"✓ MPI launchers: {', '.join(found)}"
    else:
        return False, "✗ No MPI launchers found"


def check_framework_structure() -> List[Tuple[bool, str]]:
    """Check if framework directories and files exist."""
    required_items = [
        ('core/', True),
        ('engine/', True),
        ('backends/', True),
        ('utils/', True),
        ('hpc_auto.py', True),
        ('examples/', False),
        ('tests/', False),
        ('requirements.txt', True),
    ]
    
    results = []
    base_path = Path(__file__).parent
    
    for item, required in required_items:
        path = base_path / item
        exists = path.exists()
        
        if exists:
            results.append((True, f"✓ {item}"))
        else:
            status = "✗" if required else "○"
            note = " (required)" if required else " (optional)"
            results.append((not required, f"{status} {item}{note}"))
    
    return results


def main():
    """Run all system checks."""
    print("="*70)
    print("HPC-ScaleTest System Check")
    print("="*70)
    print()
    
    all_passed = True
    
    # Python version
    print("Python Environment:")
    passed, msg = check_python_version()
    print(f"  {msg}")
    all_passed &= passed
    print()
    
    # Required Python modules
    print("Required Python Modules:")
    modules = [
        ('yaml', True),  # pyyaml
        ('pathlib', True),
        ('subprocess', True),
        ('argparse', True),
    ]
    
    for module, required in modules:
        passed, msg = check_module(module, required)
        print(f"  {msg}")
        all_passed &= passed
    print()
    
    # Optional Python modules
    print("Optional Python Modules:")
    optional_modules = [
        ('pytest', False),
        ('matplotlib', False),
        ('numpy', False),
    ]
    
    for module, required in optional_modules:
        passed, msg = check_module(module, required)
        print(f"  {msg}")
    print()
    
    # HPC Environment
    print("HPC Environment:")
    
    # Scheduler
    passed, scheduler_name, msg = check_scheduler()
    print(f"  {msg}")
    all_passed &= passed
    
    # MPI
    passed, msg = check_mpi()
    print(f"  {msg}")
    all_passed &= passed
    print()
    
    # Module system
    print("Module System:")
    module_cmds = [
        'module',
        'modulecmd',
    ]
    
    module_found = False
    for cmd in module_cmds:
        passed, msg = check_command(cmd, required=False)
        if passed:
            module_found = True
        print(f"  {msg}")
    
    if not module_found:
        print("  ○ No module system detected (optional)")
    print()
    
    # Framework structure
    print("Framework Structure:")
    structure_results = check_framework_structure()
    for passed, msg in structure_results:
        print(f"  {msg}")
        all_passed &= passed
    print()
    
    # Framework imports
    print("Framework Imports:")
    framework_modules = [
        ('core.config', True),
        ('engine.orchestrator', True),
        ('utils.config_parser', True),
        ('backends.schedulers.slurm', True),
    ]
    
    for module, required in framework_modules:
        passed, msg = check_module(module, required)
        print(f"  {msg}")
        all_passed &= passed
    print()
    
    # Summary
    print("="*70)
    if all_passed:
        print("✓ System check PASSED - All required components are available")
        print()
        print("You can now run HPC-ScaleTest:")
        print("  python3 hpc_auto.py --help")
        return 0
    else:
        print("✗ System check FAILED - Some required components are missing")
        print()
        print("Please install missing dependencies:")
        print("  pip install --user -r requirements.txt")
        return 1


if __name__ == '__main__':
    sys.exit(main())
