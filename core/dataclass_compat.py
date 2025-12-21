#!/usr/bin/env python3
"""
Dataclasses backport compatibility for Python 3.6

This module provides a minimal dataclass-like functionality for Python 3.6
by using simple class generation.
"""

import sys

# Check if dataclasses is available (Python 3.7+)
try:
    from dataclasses import dataclass, field
    HAS_DATACLASSES = True
except ImportError:
    HAS_DATACLASSES = False
    
    def field(default=None, default_factory=None):
        """Minimal field implementation for Python 3.6."""
        if default_factory is not None:
            return default_factory
        return default
    
    def dataclass(cls):
        """
        Minimal dataclass decorator for Python 3.6.
        
        This is a simplified version that just returns the class as-is.
        The class should define __init__ manually if needed.
        """
        # Add a simple __repr__ if not already defined
        if not hasattr(cls, '__repr__'):
            def __repr__(self):
                attrs = []
                for attr in dir(self):
                    if not attr.startswith('_') and not callable(getattr(self, attr)):
                        attrs.append("{}={}".format(attr, repr(getattr(self, attr))))
                return "{}({})".format(cls.__name__, ', '.join(attrs))
            cls.__repr__ = __repr__
        
        return cls

__all__ = ['dataclass', 'field', 'HAS_DATACLASSES']
