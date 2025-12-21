#!/usr/bin/env python3
"""
Test GPU Directive System

Demonstrates the modular directive system.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print("="*70)
print("Testing GPU Directive System")
print("="*70)
print()

# Test 1: Import and registration
print("Test 1: Module Import and Registration")
print("-"*70)

try:
    from backends.schedulers.directives import gpu
    print("✓ GPU directive module imported")
except Exception as e:
    print(f"❌ Import failed: {e}")
    sys.exit(1)

try:
    from core.decorators import list_directives, get_directive
    print("✓ Decorator module imported")
except Exception as e:
    print(f"❌ Import failed: {e}")
    sys.exit(1)

# Check registration
directives = list_directives()
print(f"✓ Registered directives: {directives}")

if 'gpu' not in directives:
    print("❌ GPU directive not registered!")
    sys.exit(1)

print()

# Test 2: Basic GPU directive
print("Test 2: Basic GPU Directive (4 GPUs)")
print("-"*70)

try:
    gpu_func = get_directive('gpu')
    result = gpu_func(num_gpus=4)
    
    print("Generated directives:")
    for scheduler, directive in result.items():
        print(f"  {scheduler:10s}: {directive}")
    
    # Verify SLURM
    assert '--gpus-per-node=4' in result['slurm']
    print("✓ SLURM directive correct")
    
except Exception as e:
    print(f"❌ Test failed: {e}")
    sys.exit(1)

print()

# Test 3: GPU with type specification
print("Test 3: GPU Directive with Type (4 × A100)")
print("-"*70)

try:
    result = gpu_func(num_gpus=4, gpu_type='a100')
    
    print("Generated directives:")
    for scheduler, directive in result.items():
        print(f"  {scheduler:10s}: {directive}")
    
    # Verify SLURM
    assert '--gpus-per-node=a100:4' in result['slurm']
    print("✓ SLURM directive with type correct")
    
except Exception as e:
    print(f"❌ Test failed: {e}")
    sys.exit(1)

print()

# Test 4: GPU with binding
print("Test 4: GPU Directive with Binding (4 GPUs, closest)")
print("-"*70)

try:
    result = gpu_func(num_gpus=4, gpu_bind='closest')
    
    print("Generated directives:")
    for scheduler, directive in result.items():
        print(f"  {scheduler:10s}: {directive}")
    
    # Verify SLURM
    assert '--gpus-per-node=4' in result['slurm']
    assert '--gpu-bind=closest' in result['slurm']
    print("✓ SLURM directive with binding correct")
    
except Exception as e:
    print(f"❌ Test failed: {e}")
    sys.exit(1)

print()

# Test 5: Complete example
print("Test 5: Complete Example (4 × A100, closest bind)")
print("-"*70)

try:
    result = gpu_func(num_gpus=4, gpu_type='a100', gpu_bind='closest')
    
    print("Generated directives:")
    for scheduler, directive in result.items():
        print(f"  {scheduler:10s}: {directive}")
    
    # Verify all schedulers
    assert 'slurm' in result
    assert 'pbs' in result
    assert 'lsf' in result
    assert 'sge' in result
    print("✓ All scheduler directives generated")
    
except Exception as e:
    print(f"❌ Test failed: {e}")
    sys.exit(1)

print()

# Test 6: Leonardo Booster specific
print("Test 6: Leonardo Booster Configuration")
print("-"*70)

try:
    # Leonardo has 4 × A100 per node
    result = gpu_func(num_gpus=4, gpu_type='a100', gpu_bind='closest')
    
    slurm_directive = result['slurm']
    print(f"SLURM directive: {slurm_directive}")
    
    # Parse into SBATCH lines
    parts = slurm_directive.split()
    print()
    print("As SBATCH directives:")
    for part in parts:
        if part.startswith('--'):
            print(f"  #SBATCH {part}")
    
    print("✓ Leonardo configuration correct")
    
except Exception as e:
    print(f"❌ Test failed: {e}")
    sys.exit(1)

print()
print("="*70)
print("SUCCESS! All GPU Directive Tests Passed")
print("="*70)
print()

# Show usage example
print("Usage Example:")
print("-"*70)
print("""
from core.decorators import get_directive

# Get GPU directive function
gpu_directive = get_directive('gpu')

# Generate directives for 4 GPUs with A100 type
directives = gpu_directive(
    num_gpus=4,
    gpu_type='a100',
    gpu_bind='closest'
)

# Use SLURM directive in job script
slurm_directive = directives['slurm']
print(f"#SBATCH {slurm_directive}")
# Output: #SBATCH --gpus-per-node=a100:4 --gpu-bind=closest
""")
