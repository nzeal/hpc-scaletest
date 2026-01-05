"""
Job Generator Modules for HPC-ScaleTest

This package provides modular, hardware-specific job script generators.
CPU and GPU job generation are separated for clarity and maintainability.

Design Principles:
==================
- MODULAR: Separate files for CPU and GPU job generation
- NO HARDCODED VALUES: All parameters from configuration
- SYSTEM AGNOSTIC: Works on any HPC system
- APPLICATION AGNOSTIC: Works with any MPI application

Usage:
======
    from engine.job_generators import CPUJobGenerator, GPUJobGenerator
    
    # CPU job
    cpu_gen = CPUJobGenerator(config)
    script = cpu_gen.generate_script(num_nodes=4)
    
    # GPU job
    gpu_gen = GPUJobGenerator(config)
    script = gpu_gen.generate_script(num_nodes=4)

Author: HPC-ScaleTest Contributors
"""

from .cpu_job_generator import CPUJobGenerator
from .gpu_job_generator import GPUJobGenerator
from .base import JobGeneratorBase, JobGeneratorConfig

__all__ = [
    'CPUJobGenerator',
    'GPUJobGenerator', 
    'JobGeneratorBase',
    'JobGeneratorConfig'
]
