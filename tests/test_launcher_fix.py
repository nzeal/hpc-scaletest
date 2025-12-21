#!/usr/bin/env python3
"""
Quick test for launcher instantiation.
Verifies that abstract methods are properly implemented.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print("Testing launcher instantiation...")
print()

try:
    from backends.launchers.mpirun import MpiRunLauncher
    print("✓ MpiRunLauncher imported successfully")
    
    # Try to instantiate
    launcher = MpiRunLauncher()
    print("✓ MpiRunLauncher instantiated successfully")
    
    # Check method exists
    supports = launcher.supports_gpu_binding()
    print("✓ supports_gpu_binding() works: {}".format(supports))
    
except Exception as e:
    print("❌ MpiRunLauncher failed: {}".format(e))
    sys.exit(1)

print()

try:
    from backends.launchers.srun import SrunLauncher
    print("✓ SrunLauncher imported successfully")
    
    # Try to instantiate
    launcher = SrunLauncher()
    print("✓ SrunLauncher instantiated successfully")
    
    # Check method exists
    supports = launcher.supports_gpu_binding()
    print("✓ supports_gpu_binding() works: {}".format(supports))
    
except Exception as e:
    print("❌ SrunLauncher failed: {}".format(e))
    sys.exit(1)

print()
print("="*70)
print("SUCCESS! All launchers instantiate correctly")
print("="*70)
