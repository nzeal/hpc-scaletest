#!/usr/bin/env python3
"""
Pre-Submission Configuration Validator

Validates HPC job configurations before submission to catch common issues:
- Invalid partition names
- GPU/CPU mismatch
- Incorrect MPI task counts
- Missing required parameters
"""

import logging
from typing import Dict, List, Tuple, Any
from pathlib import Path

from utils.gpu_detection import GPUDetector
from utils.partition_validator import PartitionValidator
from utils.hardware_detection import HardwareDetector


logger = logging.getLogger(__name__)


class ConfigValidator:
    """Validate job configuration before submission."""
    
    def __init__(self):
        self.gpu_detector = GPUDetector()
        self.partition_validator = PartitionValidator()
        self.hardware_detector = HardwareDetector()
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.info: List[str] = []
    
    def validate(self, config: Dict[str, Any]) -> Tuple[bool, List[str], List[str]]:
        """
        Validate complete configuration.
        
        Args:
            config: Configuration dictionary
            
        Returns:
            (is_valid, errors, warnings)
        """
        self.errors = []
        self.warnings = []
        self.info = []
        
        # Required field validation
        self._validate_required_fields(config)
        
        # Partition validation
        self._validate_partition(config)
        
        # Hardware validation
        self._validate_hardware(config)
        
        # Scaling validation
        self._validate_scaling(config)
        
        # Resource validation
        self._validate_resources(config)
        
        # Path validation
        self._validate_paths(config)
        
        is_valid = len(self.errors) == 0
        return is_valid, self.errors, self.warnings
    
    def _validate_required_fields(self, config: Dict[str, Any]):
        """Validate required configuration fields."""
        required = ['repository', 'scaling', 'nodes', 'hardware', 'partition', 'account']
        
        for field in required:
            if field not in config or not config[field]:
                self.errors.append(f"Missing required field: '{field}'")
    
    def _validate_partition(self, config: Dict[str, Any]):
        """Validate partition name and accessibility."""
        partition = config.get('partition')
        if not partition:
            return
        
        is_valid, message = self.partition_validator.validate_partition(partition)
        
        if not is_valid:
            self.errors.append(f"Partition validation failed:\n{message}")
        else:
            self.info.append(message)
            
            # Get partition info for further validation
            partition_info = self.partition_validator.get_partition_info(partition)
            if partition_info:
                # Check if partition has GPUs when GPU hardware is requested
                hardware_type = config.get('hardware', 'cpu')
                if hardware_type == 'gpu' and partition_info.gpus_per_node == 0:
                    self.warnings.append(
                        f"Partition '{partition}' has 0 GPUs but hardware type is 'gpu'. "
                        f"This may cause job failures."
                    )
    
    def _validate_hardware(self, config: Dict[str, Any]):
        """Validate hardware configuration."""
        hardware_type = config.get('hardware', 'cpu')
        partition = config.get('partition')
        
        if hardware_type not in ['cpu', 'gpu']:
            self.errors.append(f"Invalid hardware type: '{hardware_type}'. Must be 'cpu' or 'gpu'")
            return
        
        if hardware_type == 'gpu':
            # Try to detect GPU configuration
            gpu_info = self.gpu_detector.detect(partition_name=partition)
            
            if gpu_info:
                self.info.append(
                    f"✓ Detected {gpu_info.count_per_node} {gpu_info.vendor.upper()} "
                    f"GPU(s) per node ({gpu_info.model})"
                )
                
                # Check if procs_per_node is specified and matches GPU count
                procs_per_node = config.get('procs_per_node')
                if procs_per_node and procs_per_node != gpu_info.count_per_node:
                    self.warnings.append(
                        f"procs_per_node={procs_per_node} but system has "
                        f"{gpu_info.count_per_node} GPUs. Recommend using "
                        f"procs_per_node={gpu_info.count_per_node} for optimal GPU binding."
                    )
                
                # Validate scaling configuration for GPUs
                if config.get('scaling') == 'weak':
                    initial_procs = config.get('initial_procs', [])
                    if isinstance(initial_procs, list):
                        total_procs = 1
                        for p in initial_procs:
                            total_procs *= p
                        
                        if total_procs != gpu_info.count_per_node:
                            self.warnings.append(
                                f"Initial process count ({total_procs}) doesn't match GPU count "
                                f"({gpu_info.count_per_node}). For optimal weak scaling, these should match."
                            )
            else:
                self.warnings.append(
                    "Could not detect GPUs on this system. Ensure GPU partition is accessible."
                )
    
    def _validate_scaling(self, config: Dict[str, Any]):
        """Validate scaling configuration."""
        scaling = config.get('scaling')
        
        if scaling not in ['strong', 'weak']:
            self.errors.append(f"Invalid scaling type: '{scaling}'. Must be 'strong' or 'weak'")
            return
        
        # Validate weak scaling specific parameters
        if scaling == 'weak':
            required_weak = ['initial_domain', 'initial_cells']
            for field in required_weak:
                if field not in config:
                    self.warnings.append(
                        f"Weak scaling typically requires '{field}' parameter"
                    )
            
            # Validate scaling dimensions
            dimensions = config.get('scaling_dimensions', 3)
            if dimensions not in [2, 3]:
                self.errors.append(
                    f"Invalid scaling_dimensions: {dimensions}. Must be 2 or 3"
                )
        
        # Validate initial_procs format
        initial_procs = config.get('initial_procs')
        if initial_procs:
            if isinstance(initial_procs, list):
                if len(initial_procs) != 3:
                    self.errors.append(
                        f"initial_procs must have 3 elements [x,y,z], got {len(initial_procs)}"
                    )
            elif isinstance(initial_procs, str):
                parts = initial_procs.split(',')
                if len(parts) != 3:
                    self.errors.append(
                        f"initial_procs must have 3 elements 'x,y,z', got {len(parts)}"
                    )
    
    def _validate_resources(self, config: Dict[str, Any]):
        """Validate resource requirements."""
        nodes = config.get('nodes', 1)
        
        if nodes <= 0:
            self.errors.append(f"Invalid nodes: {nodes}. Must be positive integer")
        
        if nodes > 1024:
            self.warnings.append(
                f"Large node count ({nodes}) may face scheduling delays or resource limits"
            )
        
        # Validate time limit format
        time_limit = config.get('time_limit')
        if time_limit:
            if not self._is_valid_time_format(time_limit):
                self.errors.append(
                    f"Invalid time_limit format: '{time_limit}'. "
                    f"Expected format: HH:MM:SS or DD-HH:MM:SS"
                )
    
    def _validate_paths(self, config: Dict[str, Any]):
        """Validate file paths."""
        repository = config.get('repository')
        if repository:
            # Check if it's a URL or local path
            if not (repository.startswith('http') or 
                   repository.startswith('git@') or
                   Path(repository).exists()):
                self.warnings.append(
                    f"Repository path '{repository}' not found. "
                    f"Ensure path is correct or repository is accessible."
                )
        
        # Check input file
        input_file = config.get('input_file')
        if input_file:
            # Only validate if repository is local path
            if repository and not repository.startswith(('http', 'git@')):
                full_path = Path(repository) / input_file
                if not full_path.exists():
                    self.warnings.append(
                        f"Input file '{input_file}' not found in repository"
                    )
    
    def _is_valid_time_format(self, time_str: str) -> bool:
        """Check if time string is in valid SLURM format."""
        import re
        
        # HH:MM:SS format
        pattern1 = re.compile(r'^\d{1,2}:\d{2}:\d{2}$')
        # DD-HH:MM:SS format
        pattern2 = re.compile(r'^\d{1,3}-\d{1,2}:\d{2}:\d{2}$')
        # MM:SS format
        pattern3 = re.compile(r'^\d{1,2}:\d{2}$')
        
        return bool(pattern1.match(time_str) or 
                   pattern2.match(time_str) or 
                   pattern3.match(time_str))
    
    def print_validation_report(self):
        """Print formatted validation report."""
        print("\n" + "="*70)
        print("Configuration Validation Report")
        print("="*70)
        
        if self.errors:
            print(f"\n❌ ERRORS ({len(self.errors)}):")
            for i, error in enumerate(self.errors, 1):
                print(f"  {i}. {error}")
        
        if self.warnings:
            print(f"\n⚠ WARNINGS ({len(self.warnings)}):")
            for i, warning in enumerate(self.warnings, 1):
                print(f"  {i}. {warning}")
        
        if self.info:
            print(f"\n✓ INFO ({len(self.info)}):")
            for i, info in enumerate(self.info, 1):
                print(f"  {i}. {info}")
        
        print("\n" + "="*70)
        
        if not self.errors:
            print("✓ Configuration validation PASSED")
            if self.warnings:
                print("  (with warnings - review above)")
        else:
            print("❌ Configuration validation FAILED")
            print("  Please fix errors before submitting jobs")
        
        print("="*70 + "\n")


def validate_config(config: Dict[str, Any], verbose: bool = True) -> bool:
    """
    Convenience function to validate configuration.
    
    Args:
        config: Configuration dictionary
        verbose: Print validation report
        
    Returns:
        True if valid, False otherwise
    """
    validator = ConfigValidator()
    is_valid, errors, warnings = validator.validate(config)
    
    if verbose:
        validator.print_validation_report()
    
    return is_valid


if __name__ == '__main__':
    # Test configuration validator
    import sys
    import yaml
    
    logging.basicConfig(level=logging.INFO)
    
    if len(sys.argv) < 2:
        print("Usage: python3 config_validator.py <config.yaml>")
        sys.exit(1)
    
    config_file = Path(sys.argv[1])
    if not config_file.exists():
        print(f"Error: Configuration file '{config_file}' not found")
        sys.exit(1)
    
    # Load and validate configuration
    with open(config_file) as f:
        config = yaml.safe_load(f)
    
    is_valid = validate_config(config, verbose=True)
    sys.exit(0 if is_valid else 1)
