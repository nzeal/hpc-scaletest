"""
High Bandwidth Memory (HBM) Feature Detection

Detects and configures HBM on systems that support it.
HBM is typically found on:
- Intel Knights Landing (KNL)
- Some AMD EPYC processors
- Specialized HPC systems
"""

import os
import logging
from typing import Dict, Any, List

from core.decorators import register_feature
from utils.hardware.base import HardwareFeature, DetectionError

logger = logging.getLogger(__name__)


@register_feature("hbm")
class HighBandwidthMemory(HardwareFeature):
    """
    High-bandwidth memory feature.
    
    Detects HBM-capable NUMA nodes and provides configuration
    for optimal memory placement.
    """
    
    def __init__(self):
        super().__init__()
        self.capacity_gb = 0
        self.numa_nodes = []  # List of HBM NUMA node IDs
        self.ddr_nodes = []   # List of DDR NUMA node IDs
        self.total_nodes = 0
    
    def detect(self) -> bool:
        """
        Detect HBM availability on system.
        
        Checks /sys/devices/system/node for NUMA nodes and
        identifies which ones have HBM vs DDR memory.
        
        Returns:
            bool: True if HBM nodes found
        """
        try:
            node_path = '/sys/devices/system/node'
            
            if not os.path.exists(node_path):
                logger.debug("No NUMA node information available")
                self.detected = True
                self.available = False
                return False
            
            # Check for HBM NUMA nodes
            hbm_nodes = []
            ddr_nodes = []
            
            for node_dir in os.listdir(node_path):
                if not node_dir.startswith('node'):
                    continue
                
                try:
                    node_num = int(node_dir[4:])
                except ValueError:
                    continue
                
                meminfo_path = os.path.join(node_path, node_dir, 'meminfo')
                
                if os.path.exists(meminfo_path):
                    with open(meminfo_path, 'r') as f:
                        meminfo_content = f.read()
                    
                    if self._is_hbm_node(meminfo_content):
                        hbm_nodes.append(node_num)
                        # Extract capacity
                        capacity = self._extract_memory_size(meminfo_content)
                        self.capacity_gb += capacity
                    else:
                        ddr_nodes.append(node_num)
            
            self.numa_nodes = sorted(hbm_nodes)
            self.ddr_nodes = sorted(ddr_nodes)
            self.total_nodes = len(hbm_nodes) + len(ddr_nodes)
            
            if hbm_nodes:
                self.available = True
                self.detected = True
                logger.info(f"✓ Detected HBM on NUMA nodes: {self.numa_nodes}")
                logger.info(f"  HBM capacity: {self.capacity_gb:.1f} GB")
                logger.info(f"  DDR nodes: {self.ddr_nodes}")
                return True
            else:
                self.available = False
                self.detected = True
                logger.debug("No HBM nodes detected")
                return False
        
        except Exception as e:
            logger.warning(f"HBM detection failed: {e}")
            raise DetectionError(f"Failed to detect HBM: {e}")
    
    def configure(self, bind_policy='preferred', use_memkind=True, **kwargs) -> Dict[str, Any]:
        """
        Configure HBM usage.
        
        Args:
            bind_policy: Memory binding policy
                - 'preferred': Prefer HBM but allow DDR fallback
                - 'bind': Bind strictly to HBM
                - 'interleave': Interleave across HBM nodes
            use_memkind: Use memkind library for explicit allocation
            **kwargs: Additional configuration options
        
        Returns:
            dict: Configuration with env_vars, launcher_args, etc.
        """
        config = {
            'env_vars': {},
            'launcher_args': [],
            'module_loads': [],
            'init_commands': []
        }
        
        if not self.available:
            logger.warning("HBM not available, returning empty config")
            return config
        
        # Configure based on bind policy
        if bind_policy == 'preferred':
            # Prefer HBM but allow fallback to DDR
            config['env_vars']['MEMKIND_HBW_NODES'] = \
                ','.join(map(str, self.numa_nodes))
            
            if use_memkind:
                config['env_vars']['LD_PRELOAD'] = 'libmemkind.so'
                config['module_loads'].append('memkind')
        
        elif bind_policy == 'bind':
            # Strict binding to HBM nodes
            node_list = ','.join(map(str, self.numa_nodes))
            config['launcher_args'].append(f'--membind={node_list}')
            
            # Also set environment variable
            config['env_vars']['MEMKIND_HBW_NODES'] = node_list
        
        elif bind_policy == 'interleave':
            # Interleave across HBM nodes
            node_list = ','.join(map(str, self.numa_nodes))
            config['launcher_args'].append(f'--mem-bind=interleave:{node_list}')
        
        # Additional HBM-specific settings
        if kwargs.get('huge_pages', False):
            config['env_vars']['MEMKIND_HBW_HUGETLB'] = '1'
        
        # Log configuration
        logger.info(f"✓ Configured HBM with policy: {bind_policy}")
        if config['env_vars']:
            logger.debug(f"  Environment: {config['env_vars']}")
        if config['launcher_args']:
            logger.debug(f"  Launcher args: {config['launcher_args']}")
        
        self.config = config
        return config
    
    def get_info(self) -> Dict[str, Any]:
        """Get detailed HBM information."""
        info = super().get_info()
        info.update({
            'capacity_gb': self.capacity_gb,
            'hbm_nodes': self.numa_nodes,
            'ddr_nodes': self.ddr_nodes,
            'total_numa_nodes': self.total_nodes
        })
        return info
    
    def _is_hbm_node(self, meminfo_content: str) -> bool:
        """
        Heuristic to identify HBM nodes.
        
        Checks for HBM keywords in meminfo or uses heuristics like:
        - Very high bandwidth indicators
        - Specific NUMA node ranges (e.g., nodes 4+ on KNL)
        """
        # Direct HBM indicators
        if 'HBM' in meminfo_content or 'MCDRAM' in meminfo_content:
            return True
        
        # Some systems don't label HBM explicitly
        # Additional heuristics could be added here
        # For example, checking memory size or bandwidth
        
        return False
    
    def _extract_memory_size(self, meminfo_content: str) -> float:
        """Extract memory size in GB from meminfo."""
        try:
            for line in meminfo_content.split('\n'):
                if 'MemTotal:' in line:
                    # Format: "Node 0 MemTotal:       32768000 kB"
                    parts = line.split()
                    kb = int(parts[-2])
                    return kb / (1024 * 1024)  # Convert to GB
        except Exception as e:
            logger.debug(f"Failed to extract memory size: {e}")
        
        return 0.0


# Utility functions

def get_numa_node_list() -> List[int]:
    """
    Get list of all NUMA nodes on system.
    
    Returns:
        list: NUMA node IDs
    """
    nodes = []
    node_path = '/sys/devices/system/node'
    
    if os.path.exists(node_path):
        for node_dir in os.listdir(node_path):
            if node_dir.startswith('node'):
                try:
                    node_num = int(node_dir[4:])
                    nodes.append(node_num)
                except ValueError:
                    pass
    
    return sorted(nodes)


def get_memory_info(numa_node: int) -> Dict[str, Any]:
    """
    Get memory information for specific NUMA node.
    
    Args:
        numa_node: NUMA node ID
    
    Returns:
        dict: Memory information (total, free, etc.)
    """
    info = {}
    meminfo_path = f'/sys/devices/system/node/node{numa_node}/meminfo'
    
    if os.path.exists(meminfo_path):
        with open(meminfo_path, 'r') as f:
            for line in f:
                if 'MemTotal:' in line:
                    kb = int(line.split()[-2])
                    info['total_gb'] = kb / (1024 * 1024)
                elif 'MemFree:' in line:
                    kb = int(line.split()[-2])
                    info['free_gb'] = kb / (1024 * 1024)
    
    return info
