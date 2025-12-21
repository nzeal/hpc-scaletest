"""
Hardware Feature Detection Base Classes

Provides abstract base class for hardware feature detection
and configuration.
"""

import logging
from abc import ABC, abstractmethod
from typing import Dict, List, Any

logger = logging.getLogger(__name__)


class HardwareFeature(ABC):
    """
    Abstract base class for hardware features.
    
    All hardware features (HBM, NVMe, InfiniBand, etc.) should
    inherit from this class and implement the required methods.
    """
    
    def __init__(self):
        self.available = False
        self.detected = False
        self.config = {}
    
    @abstractmethod
    def detect(self) -> bool:
        """
        Detect if feature is available on current system.
        
        Returns:
            bool: True if feature is available, False otherwise
        """
        pass
    
    @abstractmethod
    def configure(self, **params) -> Dict[str, Any]:
        """
        Configure feature for use.
        
        Args:
            **params: Feature-specific configuration parameters
        
        Returns:
            dict: Configuration dictionary with keys:
                - 'env_vars': Environment variables to set
                - 'launcher_args': Arguments for MPI launcher
                - 'module_loads': Modules to load
                - 'init_commands': Initialization commands
        """
        pass
    
    def get_info(self) -> Dict[str, Any]:
        """
        Get information about detected feature.
        
        Returns:
            dict: Feature information
        """
        return {
            'name': self.__class__.__name__,
            'available': self.available,
            'detected': self.detected,
            'config': self.config
        }
    
    def is_available(self) -> bool:
        """Check if feature is available."""
        return self.available
    
    def __repr__(self):
        return f"{self.__class__.__name__}(available={self.available})"


class DetectionError(Exception):
    """Exception raised when feature detection fails."""
    pass


class ConfigurationError(Exception):
    """Exception raised when feature configuration fails."""
    pass
