"""
Hardware Feature Discovery

Automatically discovers and configures available hardware features.
"""

import logging
from typing import Dict, List, Any
from core.registry import PluginRegistry

logger = logging.getLogger(__name__)


def discover_hardware_features(verbose=True) -> Dict[str, Any]:
    """
    Automatically discover available hardware features.
    
    Iterates through all registered features and detects which
    ones are available on the current system.
    
    Args:
        verbose: Log discovered features
    
    Returns:
        dict: Mapping of feature name to feature instance
    
    Example:
        >>> features = discover_hardware_features()
        >>> if 'hbm' in features:
        >>>     config = features['hbm'].configure(bind_policy='preferred')
    """
    detected = {}
    
    # Import all feature modules to trigger registration
    _import_feature_modules()
    
    feature_names = PluginRegistry.list_features()
    
    if verbose:
        logger.info(f"Scanning for {len(feature_names)} registered hardware features...")
    
    for feature_name in feature_names:
        try:
            # Get feature class from registry
            feature_class = PluginRegistry.get_feature(feature_name)
            
            # Instantiate feature
            feature = feature_class()
            
            # Attempt detection
            if feature.detect():
                detected[feature_name] = feature
                if verbose:
                    logger.info(f"  ✓ {feature_name}: Available")
            else:
                if verbose:
                    logger.debug(f"  ✗ {feature_name}: Not available")
        
        except Exception as e:
            logger.warning(f"  ⚠ {feature_name}: Detection failed - {e}")
    
    if verbose:
        logger.info(f"✓ Discovered {len(detected)} available feature(s)")
    
    return detected


def get_feature_configs(features: Dict[str, Any], **config_overrides) -> Dict[str, Dict]:
    """Get configurations for all detected features."""
    configs = {}
    
    for feature_name, feature in features.items():
        try:
            feature_overrides = config_overrides.get(feature_name, {})
            config = feature.configure(**feature_overrides)
            configs[feature_name] = config
        except Exception as e:
            logger.warning(f"Failed to configure {feature_name}: {e}")
    
    return configs


def merge_feature_configs(configs: Dict[str, Dict]) -> Dict[str, Any]:
    """Merge multiple feature configurations into unified config."""
    merged = {
        'env_vars': {},
        'launcher_args': [],
        'module_loads': [],
        'init_commands': []
    }
    
    for feature_name, config in configs.items():
        if 'env_vars' in config:
            merged['env_vars'].update(config['env_vars'])
        if 'launcher_args' in config:
            merged['launcher_args'].extend(config['launcher_args'])
        if 'module_loads' in config:
            merged['module_loads'].extend(config['module_loads'])
        if 'init_commands' in config:
            merged['init_commands'].extend(config['init_commands'])
    
    # Remove duplicates
    merged['launcher_args'] = list(dict.fromkeys(merged['launcher_args']))
    merged['module_loads'] = list(dict.fromkeys(merged['module_loads']))
    merged['init_commands'] = list(dict.fromkeys(merged['init_commands']))
    
    return merged


def auto_configure_hardware(verbose=True, **overrides) -> Dict[str, Any]:
    """One-step hardware detection and configuration."""
    features = discover_hardware_features(verbose=verbose)
    configs = get_feature_configs(features, **overrides)
    merged = merge_feature_configs(configs)
    
    if verbose and merged['env_vars']:
        logger.info(f"✓ Hardware environment variables: {len(merged['env_vars'])}")
    
    return merged


def _import_feature_modules():
    """Import all feature modules to trigger registration."""
    try:
        from utils.hardware import hbm
    except Exception as e:
        logger.debug(f"Failed to import some feature modules: {e}")


__all__ = [
    'discover_hardware_features',
    'get_feature_configs',
    'merge_feature_configs',
    'auto_configure_hardware'
]
