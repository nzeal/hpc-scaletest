# ====================
# core/__init__.py
# ====================
"""
Core abstractions and types for HPC-ScaleTest.
"""

from .types import (
    ScalingType,
    SchedulerBackend,
    LauncherBackend,
    ModuleBackend,
    BuildBackend,
    JobStatus
)
from .test_definition import Test
from .factory import BackendFactory
from .registry import (
    register_launcher,
    get_launcher,
    list_launchers,
    has_launcher,
    JobLauncher
)

__all__ = [
    'ScalingType',
    'SchedulerBackend',
    'LauncherBackend',
    'ModuleBackend',
    'BuildBackend',
    'JobStatus',
    'Test',
    'BackendFactory',
    'register_launcher',
    'get_launcher',
    'list_launchers',
    'has_launcher',
    'JobLauncher'
]

