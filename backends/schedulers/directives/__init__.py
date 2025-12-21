"""
Scheduler Directives Module

Auto-loads all directive modules to register them.
"""

# Import all directive modules to trigger registration
from backends.schedulers.directives import gpu

__all__ = ['gpu']
