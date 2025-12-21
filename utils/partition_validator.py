#!/usr/bin/env python3
"""
Partition Validation Utility

Validates SLURM partitions and provides helpful error messages.
"""

import subprocess
import logging
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass


logger = logging.getLogger(__name__)


@dataclass
class PartitionInfo:
    """Information about a SLURM partition."""
    name: str
    state: str
    nodes: int
    cpus_per_node: int
    gpus_per_node: int
    max_time: str
    default_qos: Optional[str] = None
    allowed_qos: List[str] = None


class PartitionValidator:
    """Validate and query SLURM partitions."""
    
    def __init__(self):
        self.available_partitions: Optional[List[PartitionInfo]] = None
    
    def get_available_partitions(self) -> List[PartitionInfo]:
        """
        Query available SLURM partitions.
        
        Returns:
            List of PartitionInfo objects
        """
        if self.available_partitions is not None:
            return self.available_partitions
        
        partitions = []
        
        try:
            result = subprocess.run(
                ['sinfo', '-o', '%R,%a,%D,%c,%G,%l', '--noheader'],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0:
                for line in result.stdout.strip().split('\n'):
                    if not line.strip():
                        continue
                    
                    parts = line.split(',')
                    if len(parts) >= 5:
                        name = parts[0].strip()
                        state = parts[1].strip()
                        nodes = int(parts[2]) if parts[2].isdigit() else 0
                        cpus = int(parts[3]) if parts[3].isdigit() else 0
                        
                        # Parse GRES for GPUs
                        gpus = 0
                        gres = parts[4] if len(parts) > 4 else ''
                        if 'gpu:' in gres:
                            import re
                            gpu_match = re.search(r'gpu:(\d+)', gres)
                            if gpu_match:
                                gpus = int(gpu_match.group(1))
                        
                        max_time = parts[5] if len(parts) > 5 else 'UNLIMITED'
                        
                        partition = PartitionInfo(
                            name=name,
                            state=state,
                            nodes=nodes,
                            cpus_per_node=cpus,
                            gpus_per_node=gpus,
                            max_time=max_time
                        )
                        partitions.append(partition)
        
        except FileNotFoundError:
            logger.warning("sinfo command not found - cannot validate partitions")
        except Exception as e:
            logger.debug(f"Failed to query partitions: {e}")
        
        self.available_partitions = partitions
        return partitions
    
    def validate_partition(self, partition_name: str) -> Tuple[bool, str]:
        """
        Validate if a partition exists and is accessible.
        
        Args:
            partition_name: Name of the partition to validate
            
        Returns:
            (is_valid, message)
        """
        partitions = self.get_available_partitions()
        
        if not partitions:
            # Cannot validate - assume it's OK
            logger.debug(f"Cannot validate partition {partition_name} - assuming valid")
            return True, f"Partition {partition_name} (validation skipped)"
        
        # Check if partition exists
        partition_names = [p.name for p in partitions]
        
        if partition_name not in partition_names:
            # Provide helpful error message with suggestions
            message = f"❌ Invalid partition: '{partition_name}'\n"
            message += f"Available partitions on this system:\n"
            for p in partitions:
                gpu_info = f", {p.gpus_per_node} GPUs/node" if p.gpus_per_node > 0 else ""
                message += f"  • {p.name} ({p.state}, {p.nodes} nodes, {p.cpus_per_node} CPUs/node{gpu_info})\n"
            
            # Try to suggest similar partition names
            similar = self._find_similar_partitions(partition_name, partition_names)
            if similar:
                message += f"\nDid you mean one of these?\n"
                for name in similar:
                    message += f"  • {name}\n"
            
            return False, message
        
        # Partition exists - check if it's usable
        partition = next(p for p in partitions if p.name == partition_name)
        
        if partition.state != 'up' and partition.state != 'UP':
            return False, f"⚠ Partition '{partition_name}' is in state: {partition.state}"
        
        return True, f"✓ Partition '{partition_name}' is valid and accessible"
    
    def get_partition_info(self, partition_name: str) -> Optional[PartitionInfo]:
        """
        Get detailed information about a partition.
        
        Args:
            partition_name: Name of the partition
            
        Returns:
            PartitionInfo object or None if not found
        """
        partitions = self.get_available_partitions()
        
        for partition in partitions:
            if partition.name == partition_name:
                return partition
        
        return None
    
    def _find_similar_partitions(self, target: str, available: List[str]) -> List[str]:
        """
        Find partition names similar to the target.
        
        Uses simple string similarity (case-insensitive substring matching).
        """
        similar = []
        target_lower = target.lower()
        
        for name in available:
            name_lower = name.lower()
            
            # Check for substring matches
            if target_lower in name_lower or name_lower in target_lower:
                similar.append(name)
                continue
            
            # Check for similar words
            target_words = set(target_lower.replace('_', ' ').replace('-', ' ').split())
            name_words = set(name_lower.replace('_', ' ').replace('-', ' ').split())
            
            if target_words & name_words:  # Common words
                similar.append(name)
        
        return similar[:3]  # Return top 3 matches
    
    def print_available_partitions(self):
        """Print a formatted table of available partitions."""
        partitions = self.get_available_partitions()
        
        if not partitions:
            print("No partition information available (sinfo not accessible)")
            return
        
        print(f"\n{'='*80}")
        print("Available SLURM Partitions:")
        print(f"{'='*80}")
        print(f"{'Name':<20} {'State':<8} {'Nodes':<8} {'CPUs':<8} {'GPUs':<8} {'MaxTime':<15}")
        print(f"{'-'*80}")
        
        for p in partitions:
            gpu_str = str(p.gpus_per_node) if p.gpus_per_node > 0 else '-'
            print(f"{p.name:<20} {p.state:<8} {p.nodes:<8} {p.cpus_per_node:<8} "
                  f"{gpu_str:<8} {p.max_time:<15}")
        
        print(f"{'='*80}\n")


def validate_partition_exists(partition_name: str) -> Tuple[bool, str]:
    """
    Convenience function to validate a partition.
    
    Args:
        partition_name: Name of the partition
        
    Returns:
        (is_valid, message)
    """
    validator = PartitionValidator()
    return validator.validate_partition(partition_name)


def list_available_partitions() -> List[str]:
    """
    Get list of available partition names.
    
    Returns:
        List of partition names
    """
    validator = PartitionValidator()
    partitions = validator.get_available_partitions()
    return [p.name for p in partitions]


if __name__ == '__main__':
    # Test partition validation
    import sys
    
    logging.basicConfig(level=logging.INFO)
    
    validator = PartitionValidator()
    validator.print_available_partitions()
    
    if len(sys.argv) > 1:
        partition = sys.argv[1]
        is_valid, message = validator.validate_partition(partition)
        print(message)
        sys.exit(0 if is_valid else 1)
