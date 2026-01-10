=============================
Extending HPC-ScaleTest
=============================

This guide shows developers how to extend HPC-ScaleTest with new scheduler directives, hardware features, and backend implementations.

Overview
========

HPC-ScaleTest uses a plugin-based architecture that allows extension without modifying core code. The system supports:

* Custom scheduler directives (GPU allocation, memory constraints, etc.)
* Hardware feature detection and configuration
* New scheduler backends (beyond SLURM/PBS)
* Custom launchers and module systems

Extension Mechanisms
====================

Registry Pattern
----------------

HPC-ScaleTest uses a registry pattern for runtime discovery:

.. code-block:: python

   # In core/registry.py
   class FeatureRegistry:
       """Central registry for hardware features and directives."""
       
       _features = {}
       _directives = {}
       
       @classmethod
       def register_feature(cls, name, feature_class):
           """Register a hardware feature."""
           cls._features[name] = feature_class
       
       @classmethod
       def register_directive(cls, name, directive_func):
           """Register a scheduler directive."""
           cls._directives[name] = directive_func
       
       @classmethod
       def get_feature(cls, name):
           """Retrieve registered feature."""
           return cls._features.get(name)
       
       @classmethod
       def get_directive(cls, name):
           """Retrieve registered directive."""
           return cls._directives.get(name)

Decorator-Based Registration
-----------------------------

The system uses decorators for automatic registration:

.. code-block:: python

   # In core/decorators.py
   from core.registry import FeatureRegistry
   
   def register_feature(name):
       """Decorator to register a hardware feature."""
       def decorator(cls):
           FeatureRegistry.register_feature(name, cls)
           return cls
       return decorator
   
   def register_directive(name):
       """Decorator to register a scheduler directive."""
       def decorator(func):
           FeatureRegistry.register_directive(name, func)
           return func
       return decorator

Adding a New Scheduler Directive
=================================

Example: GPU Allocation
-----------------------

Here's how to add GPU allocation support for SLURM:

**Step 1: Define the Directive Function**

.. code-block:: python

   # In backends/schedulers/directives/gpu.py
   from core.decorators import register_directive
   
   @register_directive("gpu")
   def gpu_directive(num_gpus, gpu_type=None, gpu_bind=None):
       """
       Generate GPU allocation directive for scheduler.
       
       Parameters
       ----------
       num_gpus : int
           Number of GPUs to allocate per node
       gpu_type : str, optional
           Specific GPU type (e.g., 'a100', 'v100')
       gpu_bind : str, optional
           GPU binding strategy
       
       Returns
       -------
       dict
           Scheduler-specific directives
       """
       directives = {
           'slurm': _slurm_gpu(num_gpus, gpu_type, gpu_bind),
           'pbs': _pbs_gpu(num_gpus, gpu_type, gpu_bind),
       }
       return directives
   
   def _slurm_gpu(num_gpus, gpu_type, gpu_bind):
       """SLURM-specific GPU directive."""
       parts = [f"--gpus-per-node={num_gpus}"]
       
       if gpu_type:
           parts[0] = f"--gpus-per-node={gpu_type}:{num_gpus}"
       
       if gpu_bind:
           parts.append(f"--gpu-bind={gpu_bind}")
       
       return ' '.join(parts)

**Step 2: Use in Application**

.. code-block:: python

   from engine.job_builder import JobBuilder
   
   builder = JobBuilder(scheduler_type='slurm')
   builder.add_directive('gpu', num_gpus=4, gpu_type='a100')
   
   script = builder.build_script(
       job_name='gpu_job',
       commands=['./my_app --use-gpu']
   )

Adding a Hardware Feature
==========================

Example: High-Bandwidth Memory
-------------------------------

**Step 1: Define Feature Class**

.. code-block:: python

   # In utils/hardware/hbm.py
   from abc import ABC, abstractmethod
   from core.decorators import register_feature
   
   class HardwareFeature(ABC):
       """Base class for hardware features."""
       
       @abstractmethod
       def detect(self):
           """Detect if feature is available."""
           pass
       
       @abstractmethod
       def configure(self, **params):
           """Configure feature for use."""
           pass
   
   @register_feature("hbm")
   class HighBandwidthMemory(HardwareFeature):
       """High-bandwidth memory feature."""
       
       def __init__(self):
           self.available = False
           self.capacity_gb = 0
           self.numa_nodes = []
       
       def detect(self):
           """Detect HBM availability on system."""
           import os
           
           node_path = '/sys/devices/system/node'
           if not os.path.exists(node_path):
               return False
           
           # Check for HBM NUMA nodes
           hbm_nodes = []
           for node_dir in os.listdir(node_path):
               if not node_dir.startswith('node'):
                   continue
               
               meminfo = os.path.join(node_path, node_dir, 'meminfo')
               if os.path.exists(meminfo):
                   with open(meminfo, 'r') as f:
                       if self._is_hbm_node(f.read()):
                           node_num = int(node_dir[4:])
                           hbm_nodes.append(node_num)
           
           if hbm_nodes:
               self.available = True
               self.numa_nodes = hbm_nodes
               return True
           
           return False
       
       def configure(self, bind_policy='preferred'):
           """Configure HBM usage."""
           config = {
               'env_vars': {},
               'launcher_args': []
           }
           
           if not self.available:
               return config
           
           if bind_policy == 'preferred':
               config['env_vars']['MEMKIND_HBW_NODES'] = \
                   ','.join(map(str, self.numa_nodes))
           
           elif bind_policy == 'bind':
               node_list = ','.join(map(str, self.numa_nodes))
               config['launcher_args'].append(
                   f'--membind={node_list}'
               )
           
           return config
       
       def _is_hbm_node(self, meminfo_content):
           """Heuristic to identify HBM nodes."""
           return 'HBM' in meminfo_content

**Step 2: Runtime Discovery**

.. code-block:: python

   # In utils/hardware/__init__.py
   from core.registry import FeatureRegistry
   
   def discover_hardware_features():
       """Automatically discover available hardware features."""
       detected = {}
       
       for feature_name in FeatureRegistry.list_features():
           feature_class = FeatureRegistry.get_feature(feature_name)
           feature = feature_class()
           
           if feature.detect():
               detected[feature_name] = feature
       
       return detected

**Step 3: Use in Job Configuration**

.. code-block:: python

   from utils.hardware import discover_hardware_features
   
   # Detect hardware
   features = discover_hardware_features()
   
   if 'hbm' in features:
       hbm = features['hbm']
       config = hbm.configure(bind_policy='preferred')
       
       # Use config in job setup
       job_env = config['env_vars']
       launcher_args = config['launcher_args']

Complete Minimal Example
=========================

Here's a complete example adding a custom memory directive:

.. code-block:: python

   """
   custom_memory.py - Custom memory allocation directive
   
   Usage:
       1. Place this file in your project
       2. Import it before creating jobs
       3. Use the directive in job builder
   """
   
   from core.decorators import register_directive
   from engine.job_builder import JobBuilder
   
   @register_directive("memory_per_node")
   def memory_per_node(memory_gb):
       """
       Allocate memory per node.
       
       Parameters
       ----------
       memory_gb : int
           Memory in GB per node
       
       Returns
       -------
       dict
           Scheduler directives
       """
       return {
           'slurm': f"--mem={memory_gb}G",
           'pbs': f"-l mem={memory_gb}gb",
       }
   
   # Example usage
   if __name__ == '__main__':
       # Import the directive (registers it automatically)
       import custom_memory
       
       # Create job builder
       builder = JobBuilder(scheduler_type='slurm')
       
       # Add custom directive
       builder.add_directive('memory_per_node', memory_gb=128)
       
       # Build script
       script = builder.build_script(
           job_name='memory_test',
           commands=['./my_app']
       )
       
       print(script)
       # Output:
       # #!/bin/bash
       # #SBATCH --job-name=memory_test
       # #SBATCH --mem=128G
       # 
       # ./my_app

Testing Extensions
==================

Unit Test Example
-----------------

.. code-block:: python

   # test_custom_directive.py
   import unittest
   from core.registry import FeatureRegistry
   from custom_memory import memory_per_node
   
   class TestMemoryDirective(unittest.TestCase):
       
       def test_registration(self):
           """Test directive is registered."""
           directive = FeatureRegistry.get_directive('memory_per_node')
           self.assertIsNotNone(directive)
       
       def test_slurm_output(self):
           """Test SLURM directive format."""
           result = memory_per_node(memory_gb=64)
           self.assertEqual(result['slurm'], '--mem=64G')
       
       def test_pbs_output(self):
           """Test PBS directive format."""
           result = memory_per_node(memory_gb=64)
           self.assertEqual(result['pbs'], '-l mem=64gb')
   
   if __name__ == '__main__':
       unittest.main()

Best Practices
==============

1. **Use Decorators for Registration**
   
   Automatic registration via ``@register_feature`` and ``@register_directive``

2. **Support Multiple Schedulers**
   
   Return dictionary with scheduler-specific implementations

3. **Provide Sensible Defaults**
   
   Make parameters optional when reasonable

4. **Document Thoroughly**
   
   Include docstrings with parameter descriptions

5. **Write Tests**
   
   Verify registration and functionality

6. **Handle Errors Gracefully**
   
   Validate inputs and provide clear messages

Extension Points Summary
========================

**Scheduler Directives**
   Register via ``@register_directive``

**Hardware Features**
   Extend ``HardwareFeature`` and use ``@register_feature``

**Scheduler Backends**
   Extend ``AbstractScheduler`` and register with factory

**Runtime Discovery**
   Features auto-discovered through registry
