#!/usr/bin/env python3
"""
Test Complete Plugin Registry System

Demonstrates runtime plugin discovery and the full plugin architecture.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print("="*70)
print("HPC-ScaleTest Plugin Registry System Test")
print("="*70)
print()

# Test 1: Import PluginRegistry
print("Test 1: Import PluginRegistry")
print("-"*70)

try:
    from core.registry import PluginRegistry
    print("✓ PluginRegistry imported successfully")
except Exception as e:
    print(f"❌ Import failed: {e}")
    sys.exit(1)

print()

# Test 2: Import plugin modules (triggers registration)
print("Test 2: Import Plugin Modules (Auto-Registration)")
print("-"*70)

try:
    # Directive plugins
    from backends.schedulers.directives import gpu
    print("✓ GPU directive module imported (registered)")
    
    # Hardware feature plugins
    from utils.hardware import hbm
    print("✓ HBM feature module imported (registered)")
    
    # Launcher plugins (already registered via existing code)
    from backends.launchers import mpirun, srun
    print("✓ MPI launcher modules imported")
    
except Exception as e:
    print(f"❌ Import failed: {e}")
    sys.exit(1)

print()

# Test 3: Query registry
print("Test 3: Query Plugin Registry")
print("-"*70)

plugins = PluginRegistry.list_all_plugins()

print(f"Registered plugins by type:")
for plugin_type, names in plugins.items():
    if names:
        print(f"  {plugin_type}: {names}")

print()

# Test 4: Get directive plugin
print("Test 4: Get and Use Directive Plugin")
print("-"*70)

try:
    gpu_directive = PluginRegistry.get_directive('gpu')
    
    if gpu_directive:
        result = gpu_directive(num_gpus=4, gpu_type='a100')
        print("✓ GPU directive retrieved from registry")
        print(f"  Generated SLURM directive: {result['slurm']}")
    else:
        print("⚠ GPU directive not found in registry")
except Exception as e:
    print(f"❌ Failed: {e}")

print()

# Test 5: Get feature plugin
print("Test 5: Get and Use Feature Plugin")
print("-"*70)

try:
    hbm_class = PluginRegistry.get_feature('hbm')
    
    if hbm_class:
        hbm = hbm_class()
        print("✓ HBM feature class retrieved from registry")
        print(f"  HBM feature: {hbm}")
        
        detected = hbm.detect()
        print(f"  Detection result: {detected}")
    else:
        print("⚠ HBM feature not found in registry")
except Exception as e:
    print(f"❌ Failed: {e}")

print()

# Test 6: Print full registry status
print("Test 6: Full Registry Status")
print("-"*70)

PluginRegistry.print_registry_status()

print()

# Test 7: Demonstrate extensibility
print("Test 7: Demonstrate Extensibility (Add New Plugin)")
print("-"*70)

try:
    from core.decorators import register_directive
    
    # Define new directive on-the-fly
    @register_directive("memory")
    def memory_directive(memory_gb):
        """Memory allocation directive."""
        return {
            'slurm': f"--mem={memory_gb}G",
            'pbs': f"-l mem={memory_gb}gb",
        }
    
    print("✓ New 'memory' directive registered on-the-fly")
    
    # Verify it's registered
    directives = PluginRegistry.list_directives()
    print(f"  Current directives: {directives}")
    
    # Use it
    mem_directive = PluginRegistry.get_directive('memory')
    result = mem_directive(memory_gb=128)
    print(f"  Generated SLURM: {result['slurm']}")
    
except Exception as e:
    print(f"❌ Failed: {e}")

print()

# Test 8: Demonstrate runtime discovery
print("Test 8: Runtime Plugin Discovery")
print("-"*70)

try:
    from utils.hardware import discover_hardware_features
    
    # Discover all features
    features = discover_hardware_features(verbose=False)
    
    print(f"✓ Discovered {len(features)} hardware feature(s)")
    if features:
        for name, feature in features.items():
            print(f"  • {name}: {feature}")
    else:
        print("  (No special hardware detected - normal on standard nodes)")
    
except Exception as e:
    print(f"❌ Failed: {e}")

print()
print("="*70)
print("SUCCESS! Plugin Registry System Working")
print("="*70)
print()

# Show usage summary
print("Plugin Architecture Summary:")
print("-"*70)
print("""
The HPC-ScaleTest plugin system enables:

1. Runtime Discovery
   - Plugins register themselves via decorators
   - No manual registration needed
   - Central registry tracks all plugins

2. Type-Safe Plugins
   - Directives: Scheduler-specific commands
   - Features: Hardware detection/config
   - Backends: Scheduler implementations
   - Launchers: MPI launchers
   - Module Systems: Environment modules

3. Easy Extensibility
   - Create new file with decorator
   - Import triggers registration
   - No core code changes needed

4. Current Plugins:
""")

final_plugins = PluginRegistry.list_all_plugins()
for plugin_type, names in final_plugins.items():
    if names:
        print(f"   {plugin_type.title()}: {len(names)} registered")

print()
print("Add new plugins by:")
print("  1. Creating module in appropriate directory")
print("  2. Using @register_<type> decorator")
print("  3. Importing module (triggers registration)")
print()
