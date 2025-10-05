"""
__init__.py files for all packages
"""

# core/__init__.py
"""Core components for HPC-ScaleTest."""

from .abstracts import (
    AbstractScheduler,
    AbstractLauncher,
    AbstractModuleSystem,
    AbstractBuildSystem,
    LaunchResult,
    JobResults
)
from .types import (
    ScalingType,
    JobStatus,
    BackendType,
    SchedulerBackend,
    LauncherBackend,
    ModuleBackend,
    BuildBackend,
    JobID,
    BuildPath,
    InstallPath
)
from .config import (
    BackendConfig,
    ResourceConfig,
    ScalingConfig,
    TestConfig
)
from .test_definition import Test, test
from .factory import BackendFactory

__all__ = [
    'AbstractScheduler',
    'AbstractLauncher',
    'AbstractModuleSystem',
    'AbstractBuildSystem',
    'LaunchResult',
    'JobResults',
    'ScalingType',
    'JobStatus',
    'BackendType',
    'SchedulerBackend',
    'LauncherBackend',
    'ModuleBackend',
    'BuildBackend',
    'JobID',
    'BuildPath',
    'InstallPath',
    'BackendConfig',
    'ResourceConfig',
    'ScalingConfig',
    'TestConfig',
    'Test',
    'test',
    'BackendFactory'
]


# backends/__init__.py
"""Backend implementations for HPC-ScaleTest."""

__all__ = []


# backends/schedulers/__init__.py
"""Scheduler backend implementations."""

from .local import LocalScheduler
from .slurm import SlurmScheduler

__all__ = ['LocalScheduler', 'SlurmScheduler']


# backends/launchers/__init__.py
"""Launcher backend implementations."""

from .srun import SrunLauncher
from .mpirun import MpiRunLauncher

__all__ = ['SrunLauncher', 'MpiRunLauncher']


# backends/modules/__init__.py
"""Module system backend implementations."""

from .nomod import NoModBackend
from .tmod import TModBackend
from .tmod4 import TMod4Backend
from .lmod import LModBackend

__all__ = ['NoModBackend', 'TModBackend', 'TMod4Backend', 'LModBackend']


# backends/builds/__init__.py
"""Build system backend implementations."""

from .make import MakeBackend
from .cmake import CMakeBackend
from .autotools import AutotoolsBackend
from .easybuild import EasyBuildBackend
from .spack import SpackBackend

__all__ = [
    'MakeBackend',
    'CMakeBackend',
    'AutotoolsBackend',
    'EasyBuildBackend',
    'SpackBackend'
]


# engine/__init__.py
"""Test execution engine."""

from .scaling import ScalingEngine, NodeConfig
from .runner import TestRunner

__all__ = ['ScalingEngine', 'NodeConfig', 'TestRunner']


# utils/__init__.py
"""Utility functions."""

from .file_utils import (
    modify_file,
    copy_file,
    create_directory,
    write_file,
    read_file
)
from .logging_config import setup_logging

__all__ = [
    'modify_file',
    'copy_file',
    'create_directory',
    'write_file',
    'read_file',
    'setup_logging'
]


# tests/__init__.py
"""Unit tests for HPC-ScaleTest."""

__all__ = []
