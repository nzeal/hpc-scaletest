"""
Dynamic system resource detection for HPC environments.

This module automatically detects:
- Total nodes available
- CPU cores per node
- GPUs per node
- Memory per node
- Scheduler type
"""

import subprocess
import logging
import re
import os
from typing import Dict, Optional, Tuple, List, Any
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class SystemResources:
    """Detected system resources."""
    total_nodes: int = 1
    cores_per_node: int = 1
    gpus_per_node: int = 0
    memory_per_node_gb: float = 0.0
    scheduler_type: str = "local"
    partition_name: Optional[str] = None
    node_list: List[str] = field(default_factory=list)


class SystemDetector:
    """Detects HPC system resources dynamically."""
    
    def __init__(self):
        """Initialize system detector."""
        self.resources = SystemResources()
    
    def detect_all(self, partition: Optional[str] = None) -> SystemResources:
        """
        Detect all system resources.
        
        Args:
            partition: Optional partition name to query specific partition
            
        Returns:
            SystemResources object with detected values
        """
        logger.info("Detecting system resources...")
        
        # Detect scheduler
        self.resources.scheduler_type = self._detect_scheduler()
        logger.info(f"  Scheduler: {self.resources.scheduler_type}")
        
        # Detect based on scheduler type
        if self.resources.scheduler_type == "slurm":
            self._detect_slurm_resources(partition)
        elif self.resources.scheduler_type == "pbs":
            self._detect_pbs_resources(partition)
        else:
            # Fall back to local system detection
            self._detect_local_resources()
        
        logger.info(f"  Total nodes: {self.resources.total_nodes}")
        logger.info(f"  Cores per node: {self.resources.cores_per_node}")
        logger.info(f"  GPUs per node: {self.resources.gpus_per_node}")
        logger.info(f"  Memory per node: {self.resources.memory_per_node_gb:.1f} GB")
        
        return self.resources
    
    def _detect_scheduler(self) -> str:
        """
        Detect which scheduler is available.
        
        Returns:
            Scheduler type: 'slurm', 'pbs', or 'local'
        """
        # Check for Slurm
        if self._command_exists('sinfo'):
            return 'slurm'
        
        # Check for PBS/Torque
        if self._command_exists('qstat'):
            return 'pbs'
        
        # Default to local
        return 'local'
    
    def _detect_slurm_resources(self, partition: Optional[str] = None):
        """
        Detect resources on Slurm systems.
        
        Args:
            partition: Optional partition to query
        """
        # Use new SlurmDetector when available for robust detection
        try:
            from utils.slurm_detector import get_slurm_detector

            detector = get_slurm_detector()

            # If user supplied a partition, try to use it; otherwise auto-select
            if partition:
                selected = partition
            else:
                # Default to DCGP if it fits, otherwise auto-select
                selected = detector.select_best_partition(self.resources.total_nodes or 1,
                                                          self.resources.cores_per_node or detector.get_node_cores(),
                                                          120) or next(iter(detector.partitions.keys()), None)

            # Fill resources from detector
            self.resources.scheduler_type = 'slurm'
            self.resources.partition_name = selected
            self.resources.total_nodes = sum(p.total_nodes for p in detector.partitions.values())
            self.resources.cores_per_node = int(detector.get_node_cores(selected) or 112)
            self.resources.memory_per_node_gb = float(detector.get_node_memory_gb(selected) or 502.0)
            # GPU detection
            # pick max gpu_per_node across partitions as a safe default
            self.resources.gpus_per_node = max(p.gpu_per_node for p in detector.partitions.values() if p.gpu_per_node >= 0)
            # Collect node list (names) from detector (sample up to 10)
            self.resources.node_list = list(detector.nodes.keys())[:50]
            logger.info(f"Auto-selected partition: {self.resources.partition_name}")
            return
        except Exception as e:
            logger.warning(f"Slurm detector unavailable or failed: {e}")

        # Fallback behavior if SlurmDetector not available
        if partition:
            # Try to find exact or partial match for partition
            actual_partition = self._resolve_partition_name(partition)
            if actual_partition != partition:
                logger.info(f"Resolved '{partition}' to '{actual_partition}'")
            partition = actual_partition
            self.resources.partition_name = partition
            logger.info(f"Querying Slurm partition: {partition}")
        
        # Try multiple detection methods in order of reliability
        success = False
        
        # Method 1: sinfo with GRES (most common)
        logger.info("Trying Method 1: sinfo with GRES field...")
        if self._query_sinfo_direct(partition):
            success = True
            logger.info("✓ Method 1 succeeded")
        
        # Method 2: scontrol show partition (more detailed)
        if not success or self.resources.gpus_per_node == 0:
            logger.info("Trying Method 2: scontrol show partition...")
            if self._query_scontrol_partition(partition):
                success = True
                logger.info("✓ Method 2 succeeded")
        
        # Method 3: Query individual node from partition
        if not success or self.resources.gpus_per_node == 0:
            logger.info("Trying Method 3: Node-level query...")
            if self._query_partition_node(partition):
                success = True
                logger.info("✓ Method 3 succeeded")
        
        if not success:
            raise RuntimeError(
                f"Failed to query Slurm partition '{partition}' with all methods. "
                "Ensure you are on an HPC login node with Slurm access."
            )
    
    def _resolve_partition_name(self, partition: str) -> str:
        """
        Resolve partition name - try exact match first, then partial match.
        
        Args:
            partition: Partition name (can be partial like 'dcgp' or 'booster')
            
        Returns:
            Resolved partition name
        """
        try:
            # Get list of available partitions
            cmd = ['sinfo', '-h', '-o', '%P']
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0 and result.stdout.strip():
                partitions = result.stdout.strip().split()
                # Remove asterisks (default partition marker)
                partitions = [p.replace('*', '') for p in partitions]
                partitions = list(set(partitions))  # Remove duplicates
                
                # Try exact match first
                if partition in partitions:
                    return partition
                
                # Try partial match (case-insensitive)
                partition_lower = partition.lower()
                matches = [p for p in partitions if partition_lower in p.lower()]
                
                if len(matches) == 1:
                    logger.info(f"Found matching partition: {matches[0]}")
                    return matches[0]
                elif len(matches) > 1:
                    # Prefer _usr_prod partitions
                    usr_prod_matches = [m for m in matches if '_usr_prod' in m]
                    if len(usr_prod_matches) == 1:
                        logger.info(f"Found matching partition: {usr_prod_matches[0]}")
                        return usr_prod_matches[0]
                    elif usr_prod_matches:
                        logger.info(f"Multiple matches found: {matches}. Using: {usr_prod_matches[0]}")
                        return usr_prod_matches[0]
                    else:
                        logger.info(f"Multiple matches found: {matches}. Using: {matches[0]}")
                        return matches[0]
                
                # No match - return original and let it fail with suggestions
                return partition
            
        except Exception as e:
            logger.debug(f"Failed to resolve partition name: {e}")
        
        return partition
    
    def _query_sinfo_direct(self, partition: Optional[str] = None) -> bool:
        """
        Query partition info using sinfo (direct method for HPC systems).
        
        Args:
            partition: Partition name
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Use the exact command that works: sinfo -p <partition> -o "%D %c %m %G" -h
            cmd = ['sinfo', '-o', '%D %c %m %G', '-h']
            if partition:
                cmd.insert(1, '-p')
                cmd.insert(2, partition)
            
            logger.info(f"Executing: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            logger.info(f"Return code: {result.returncode}")
            logger.info(f"stdout length: {len(result.stdout) if result.stdout else 0} chars")
            logger.info(f"stdout repr: {repr(result.stdout)}")
            
            if result.stderr:
                logger.warning(f"stderr: {result.stderr.strip()}")
            
            if result.returncode != 0:
                logger.error(f"sinfo command failed with return code {result.returncode}")
                self._suggest_available_partitions()
                return False
            
            if not result.stdout or not result.stdout.strip():
                logger.error(f"sinfo returned empty output for partition '{partition}'")
                logger.error("The partition name may be incorrect or inaccessible.")
                self._suggest_available_partitions()
                return False
            
            output = result.stdout.strip()
            logger.info(f"Parsing output: '{output}'")
            
            # Parse the output: <nodes> <cpus> <memory_MB> <GRES>
            parts = output.split()
            logger.info(f"Split into {len(parts)} parts: {parts}")
            
            if len(parts) < 3:
                logger.error(f"Insufficient fields. Expected at least 3, got {len(parts)}: {parts}")
                return False
            
            # Parse values
            try:
                self.resources.total_nodes = int(parts[0])
                self.resources.cores_per_node = int(parts[1])
                self.resources.memory_per_node_gb = float(parts[2]) / 1024.0  # MB to GB
                
                logger.info(f"✓ Nodes: {self.resources.total_nodes}")
                logger.info(f"✓ Cores/node: {self.resources.cores_per_node}")
                logger.info(f"✓ Memory/node: {self.resources.memory_per_node_gb:.1f} GB")
                
            except ValueError as e:
                logger.error(f"Failed to parse numeric values: {e}")
                return False
            
            # Parse GPU info if available
            if len(parts) >= 4:
                gpus_info = parts[3]
                logger.info(f"GPU field from sinfo: '{gpus_info}'")
                
                if gpus_info and gpus_info != '(null)' and 'gpu' in gpus_info.lower():
                    self.resources.gpus_per_node = self._parse_slurm_gres(gpus_info)
                    logger.info(f"✓ GPUs/node: {self.resources.gpus_per_node}")
                else:
                    self.resources.gpus_per_node = 0
                    if gpus_info == '(null)':
                        logger.info(f"✓ GPUs/node: 0 (GRES field is '(null)')")
                        logger.info(f"  If this partition has GPUs, they may not be configured in Slurm GRES")
                        logger.info(f"  Specify gpus_per_node explicitly in your config YAML")
                    else:
                        logger.info(f"✓ GPUs/node: 0 (no GPU info in GRES field: '{gpus_info}')")
            else:
                self.resources.gpus_per_node = 0
                logger.info(f"✓ GPUs/node: 0 (no GRES field in sinfo output)")
                logger.info(f"  sinfo only returned {len(parts)} fields, expected at least 4")
                logger.info(f"  If this partition has GPUs, specify gpus_per_node in your config")
            
            return True
            
        except subprocess.TimeoutExpired:
            logger.error("sinfo command timed out")
            return False
        except Exception as e:
            logger.error(f"sinfo query failed: {e}")
            import traceback
            logger.debug(traceback.format_exc())
            return False
    
    def _suggest_available_partitions(self):
        """List available partitions to help user."""
        try:
            logger.info("Attempting to list available partitions...")
            cmd = ['sinfo', '-h', '-o', '%P']
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0 and result.stdout.strip():
                partitions = result.stdout.strip().split()
                # Remove asterisks (default partition marker)
                partitions = [p.replace('*', '') for p in partitions]
                partitions = list(set(partitions))  # Remove duplicates
                
                logger.error(f"Available partitions: {', '.join(partitions)}")
                logger.error(f"Please use one of these partition names in your run.yaml")
            else:
                logger.error("Could not list available partitions")
        except Exception as e:
            logger.debug(f"Failed to list partitions: {e}")
    
    def _query_scontrol_partition(self, partition: Optional[str] = None) -> bool:
        """
        Query partition info using scontrol show partition.
        This method can extract GPU info from TRESPerNode.
        
        Args:
            partition: Partition name
            
        Returns:
            True if successful and got useful info, False otherwise
        """
        if not partition:
            return False
            
        try:
            cmd = ['scontrol', 'show', 'partition', partition]
            logger.info(f"Executing: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode != 0:
                logger.warning(f"scontrol command failed")
                return False
            
            if not result.stdout.strip():
                return False
            
            output = result.stdout.strip()
            logger.debug(f"scontrol output length: {len(output)} chars")
            
            # Extract TotalNodes
            nodes_match = re.search(r'TotalNodes=(\d+)', output)
            if nodes_match and self.resources.total_nodes == 0:
                self.resources.total_nodes = int(nodes_match.group(1))
                logger.info(f"  From scontrol - TotalNodes: {self.resources.total_nodes}")
            
            # Try to extract from TRESPerNode if present
            tres_match = re.search(r'TRESPerNode=([^\s]+)', output)
            if tres_match:
                tres_str = tres_match.group(1)
                logger.info(f"  Found TRESPerNode: {tres_str}")
                
                # Parse TRES format: cpu=32,mem=512000M,gres/gpu=4
                cpu_match = re.search(r'cpu=(\d+)', tres_str)
                if cpu_match:
                    self.resources.cores_per_node = int(cpu_match.group(1))
                    logger.info(f"  From TRES - CPUs/node: {self.resources.cores_per_node}")
                
                mem_match = re.search(r'mem=(\d+)([KMG]?)', tres_str)
                if mem_match:
                    mem_val = int(mem_match.group(1))
                    mem_unit = mem_match.group(2) if mem_match.group(2) else 'M'
                    if mem_unit == 'K':
                        mem_val = mem_val / 1024 / 1024
                    elif mem_unit == 'M':
                        mem_val = mem_val / 1024
                    self.resources.memory_per_node_gb = float(mem_val)
                    logger.info(f"  From TRES - Memory/node: {self.resources.memory_per_node_gb:.1f} GB")
                
                gpu_match = re.search(r'gres/gpu[:\w]*[=:](\d+)', tres_str)
                if gpu_match:
                    self.resources.gpus_per_node = int(gpu_match.group(1))
                    logger.info(f"  From TRES - GPUs/node: {self.resources.gpus_per_node}")
                    return True
            
            # If we got something useful, return True
            return self.resources.cores_per_node > 0 or self.resources.total_nodes > 0
            
        except Exception as e:
            logger.debug(f"scontrol query failed: {e}")
            return False
    
    def _query_partition_node(self, partition: Optional[str] = None) -> bool:
        """
        Query a sample node from the partition to get per-node specs.
        
        Args:
            partition: Partition name
            
        Returns:
            True if successful, False otherwise
        """
        if not partition:
            return False
        
        try:
            # Get list of nodes in this partition
            cmd = ['sinfo', '-p', partition, '-N', '-h', '-o', '%N']
            logger.info(f"Getting nodes: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode != 0 or not result.stdout.strip():
                return False
            
            nodes = [n.strip() for n in result.stdout.strip().split('\n') if n.strip()]
            if not nodes:
                return False
            
            first_node = nodes[0]
            logger.info(f"  Querying sample node: {first_node}")
            
            # Query node directly
            cmd = ['sinfo', '-n', first_node, '-o', '%c %m %G', '-h']
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode != 0 or not result.stdout.strip():
                return False
            
            output = result.stdout.strip()
            parts = output.split()
            logger.info(f"  Node query returned: {output}")
            
            if len(parts) >= 2:
                self.resources.cores_per_node = int(parts[0])
                self.resources.memory_per_node_gb = float(parts[1]) / 1024.0
                logger.info(f"  From node - Cores: {self.resources.cores_per_node}")
                logger.info(f"  From node - Memory: {self.resources.memory_per_node_gb:.1f} GB")
                
                # Parse GPU info if present
                if len(parts) >= 3:
                    gres_str = parts[2]
                    if gres_str and gres_str != '(null)':
                        self.resources.gpus_per_node = self._parse_slurm_gres(gres_str)
                        logger.info(f"  From node - GPUs: {self.resources.gpus_per_node}")
                
                # Count total nodes if needed
                if self.resources.total_nodes == 0:
                    self.resources.total_nodes = len(nodes)
                    logger.info(f"  Total nodes: {self.resources.total_nodes}")
                
                return True
            
            return False
            
        except Exception as e:
            logger.debug(f"Node query failed: {e}")
            return False

    
    def _parse_slurm_gres(self, gres_str: str) -> int:
        """
        Parse Slurm GRES string to extract GPU count.
        
        Handles various GRES formats from different systems:
        - 'gpu:4' (standard)
        - 'gpu:tesla:2' (with GPU type)
        - 'gpu:a100:4(S:0-1)' (with socket info)
        - 'gpu' (assume 1)
        - '(null)' or empty (0 GPUs)
        
        Args:
            gres_str: GRES string from sinfo
            
        Returns:
            Number of GPUs
        """
        if not gres_str or gres_str == '(null)' or gres_str == 'N/A':
            return 0
        
        # Try multiple regex patterns for different GRES formats
        patterns = [
            r'gpu[:\w]*:(\d+)',           # gpu:4 or gpu:tesla:2
            r'gpu.*\(IDX:[\d,-]+\)',      # Count GPUs from index list
            r'gpu',                        # Just 'gpu' without count
        ]
        
        for pattern in patterns:
            if pattern == r'gpu':
                # Just 'gpu' - assume 1 GPU
                if gres_str.lower() == 'gpu':
                    logger.debug(f"  GRES '{gres_str}' → assuming 1 GPU")
                    return 1
            else:
                match = re.search(pattern, gres_str, re.IGNORECASE)
                if match:
                    count = int(match.group(1))
                    logger.debug(f"  GRES '{gres_str}' → {count} GPUs (pattern: {pattern})")
                    return count
        
        # If we get here, we couldn't parse it
        logger.warning(f"  Could not parse GRES string: '{gres_str}' - assuming 0 GPUs")
        logger.warning(f"  If this system has GPUs, please specify gpus_per_node in your config")
        return 0
    
    def _detect_pbs_resources(self, partition: Optional[str] = None):
        """
        Detect resources on PBS/Torque systems.
        
        Args:
            partition: Optional queue to query
        """
        try:
            # Try pbsnodes command
            cmd = ['pbsnodes', '-a']
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0 and result.stdout.strip():
                nodes = result.stdout.split('\n\n')
                self.resources.total_nodes = len([n for n in nodes if n.strip()])
                
                # Parse first node for resources
                if nodes:
                    first_node = nodes[0]
                    
                    # Extract CPUs
                    cpu_match = re.search(r'np\s*=\s*(\d+)', first_node)
                    if cpu_match:
                        self.resources.cores_per_node = int(cpu_match.group(1))
                    
                    # Extract memory
                    mem_match = re.search(r'totmem\s*=\s*(\d+)([kmg]b)', first_node, re.IGNORECASE)
                    if mem_match:
                        mem_value = int(mem_match.group(1))
                        mem_unit = mem_match.group(2).lower()
                        if mem_unit == 'kb':
                            self.resources.memory_per_node_gb = mem_value / (1024 * 1024)
                        elif mem_unit == 'mb':
                            self.resources.memory_per_node_gb = mem_value / 1024
                        else:  # gb
                            self.resources.memory_per_node_gb = float(mem_value)
                    
                    # Extract GPUs
                    gpu_match = re.search(r'gpus\s*=\s*(\d+)', first_node)
                    if gpu_match:
                        self.resources.gpus_per_node = int(gpu_match.group(1))
            else:
                self._detect_local_resources()
                
        except Exception as e:
            logger.warning(f"PBS detection failed: {e}, using local fallback")
            self._detect_local_resources()
    
    def _detect_local_resources(self):
        """Detect resources on local system (login node)."""
        logger.info("Detecting local system resources...")
        
        # Use /proc or system commands
        self.resources.total_nodes = 1
        
        # Detect CPU cores
        try:
            # Try lscpu first
            result = subprocess.run(['lscpu'], capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                for line in result.stdout.split('\n'):
                    if line.startswith('CPU(s):'):
                        cores = line.split(':')[1].strip()
                        self.resources.cores_per_node = int(cores)
                        break
            else:
                # Fall back to /proc/cpuinfo
                with open('/proc/cpuinfo', 'r') as f:
                    self.resources.cores_per_node = sum(1 for line in f if line.startswith('processor'))
        except Exception as e:
            logger.warning(f"CPU detection failed: {e}, using default")
            self.resources.cores_per_node = os.cpu_count() or 1
        
        # Detect memory
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemTotal:'):
                        mem_kb = int(line.split()[1])
                        self.resources.memory_per_node_gb = mem_kb / (1024 * 1024)
                        break
        except Exception as e:
            logger.warning(f"Memory detection failed: {e}")
            self.resources.memory_per_node_gb = 0.0
        
        # Detect GPUs
        try:
            # Try nvidia-smi
            result = subprocess.run(['nvidia-smi', '--list-gpus'], 
                                    capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                gpu_lines = [line for line in result.stdout.split('\n') if line.strip()]
                self.resources.gpus_per_node = len(gpu_lines)
        except Exception:
            # No GPUs or nvidia-smi not available
            self.resources.gpus_per_node = 0
    
    def _command_exists(self, command: str) -> bool:
        """
        Check if a command exists in PATH.
        
        Args:
            command: Command name to check
            
        Returns:
            True if command exists
        """
        try:
            result = subprocess.run(['which', command], 
                                    capture_output=True, timeout=2)
            return result.returncode == 0
        except Exception:
            return False
    
    def get_recommended_nodes(self, max_nodes: Optional[int] = None) -> int:
        """
        Get recommended number of nodes for testing.
        
        Args:
            max_nodes: User-specified maximum nodes (optional)
            
        Returns:
            Recommended node count
        """
        if max_nodes is not None:
            # User specified, but cap at available nodes
            return min(max_nodes, self.resources.total_nodes)
        
        # Auto-select reasonable fraction of cluster
        # For small clusters, use all nodes
        if self.resources.total_nodes <= 8:
            return self.resources.total_nodes
        
        # For larger clusters, use up to 25% or 64 nodes (whichever is smaller)
        recommended = min(self.resources.total_nodes // 4, 64)
        return max(recommended, 1)
    
    def get_scaling_sequence(self, max_nodes: int) -> list:
        """
        Generate recommended node sequence for scaling tests.
        
        Args:
            max_nodes: Maximum number of nodes
            
        Returns:
            List of node counts
        """
        sequence = []
        n = 1
        
        # Powers of 2 up to max_nodes
        while n <= max_nodes:
            sequence.append(n)
            n *= 2
        
        # Add max_nodes if not already included
        if sequence[-1] != max_nodes:
            sequence.append(max_nodes)
        
        return sequence


def detect_system_resources(partition: Optional[str] = None) -> SystemResources:
    """
    Convenience function to detect system resources.
    
    Args:
        partition: Optional partition name to query
        
    Returns:
        SystemResources object
    """
    detector = SystemDetector()
    return detector.detect_all(partition)


def list_partitions() -> List[str]:
    """
    List all available Slurm partitions.
    
    Returns:
        List of partition names
    """
    try:
        cmd = ['sinfo', '-h', '-o', '%P']
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        
        if result.returncode == 0 and result.stdout.strip():
            output = result.stdout.strip().split()
            # Remove duplicates and trailing asterisks (default partition marker)
            return list({p.replace('*', '') for p in output})
    except Exception as e:
        logger.debug(f"Could not list partitions: {e}")
    
    return []


def get_partition_info(partition: str) -> Dict[str, Any]:
    """
    Get detailed hardware information for a specific partition.
    
    Args:
        partition: Partition name (e.g., 'dcgp', 'booster')
        
    Returns:
        Dictionary with partition hardware details
    """
    logger.info(f"Querying hardware info for partition: {partition}")
    
    resources = detect_system_resources(partition)
    
    info = {
        'partition_name': partition,
        'total_nodes': resources.total_nodes,
        'cores_per_node': resources.cores_per_node,
        'gpus_per_node': resources.gpus_per_node,
        'memory_per_node_gb': resources.memory_per_node_gb,
        'scheduler': resources.scheduler_type,
        'partition_type': 'GPU' if resources.gpus_per_node > 0 else 'CPU',
        'node_list': resources.node_list
    }
    
    return info


def display_partition_info(partition: str):
    """
    Display formatted hardware information for a partition.
    
    Args:
        partition: Partition name to query
    """
    info = get_partition_info(partition)
    
    print("\n" + "="*70)
    print(f"HARDWARE SPECIFICATIONS: {partition.upper()} Partition")
    print("="*70)
    print(f"Partition Type:     {info['partition_type']}")
    print(f"Total Nodes:        {info['total_nodes']}")
    print(f"CPUs per Node:      {info['cores_per_node']} cores")
    
    if info['gpus_per_node'] > 0:
        print(f"GPUs per Node:      {info['gpus_per_node']}")
    else:
        print(f"GPUs per Node:      None (CPU-only partition)")
    
    print(f"Memory per Node:    {info['memory_per_node_gb']:.1f} GB")
    print(f"Scheduler:          {info['scheduler']}")
    print("="*70)
    
    # Calculate totals
    total_cores = info['total_nodes'] * info['cores_per_node']
    total_gpus = info['total_nodes'] * info['gpus_per_node']
    total_memory = info['total_nodes'] * info['memory_per_node_gb']
    
    print(f"\nTotal Resources:")
    print(f"  Total CPU Cores:  {total_cores:,}")
    if total_gpus > 0:
        print(f"  Total GPUs:       {total_gpus:,}")
    print(f"  Total Memory:     {total_memory:,.1f} GB")
    print("="*70 + "\n")


def compare_partitions(partitions: List[str]) -> Dict[str, Dict[str, Any]]:
    """
    Compare hardware specifications across multiple partitions.
    
    Args:
        partitions: List of partition names to compare
        
    Returns:
        Dictionary mapping partition names to their hardware info
    """
    comparison = {}
    
    for partition in partitions:
        try:
            comparison[partition] = get_partition_info(partition)
        except Exception as e:
            logger.warning(f"Could not get info for partition {partition}: {e}")
            comparison[partition] = None
    
    return comparison


def display_partition_comparison(partitions: List[str]):
    """
    Display side-by-side comparison of multiple partitions.
    
    Args:
        partitions: List of partition names to compare
    """
    comparison = compare_partitions(partitions)
    
    print("\n" + "="*90)
    print("PARTITION COMPARISON")
    print("="*90)
    
    # Header
    header = f"{'Specification':<25}"
    for partition in partitions:
        header += f"{partition.upper():<20}"
    print(header)
    print("-" * 90)
    
    # Partition type
    row = f"{'Type':<25}"
    for partition in partitions:
        if comparison[partition]:
            row += f"{comparison[partition]['partition_type']:<20}"
        else:
            row += f"{'N/A':<20}"
    print(row)
    
    # Nodes
    row = f"{'Total Nodes':<25}"
    for partition in partitions:
        if comparison[partition]:
            row += f"{comparison[partition]['total_nodes']:<20}"
        else:
            row += f"{'N/A':<20}"
    print(row)
    
    # CPUs per node
    row = f"{'CPUs/Node':<25}"
    for partition in partitions:
        if comparison[partition]:
            row += f"{comparison[partition]['cores_per_node']:<20}"
        else:
            row += f"{'N/A':<20}"
    print(row)
    
    # GPUs per node
    row = f"{'GPUs/Node':<25}"
    for partition in partitions:
        if comparison[partition]:
            gpus = comparison[partition]['gpus_per_node']
            row += f"{gpus if gpus > 0 else 'None':<20}"
        else:
            row += f"{'N/A':<20}"
    print(row)
    
    # Memory per node
    row = f"{'Memory/Node (GB)':<25}"
    for partition in partitions:
        if comparison[partition]:
            row += f"{comparison[partition]['memory_per_node_gb']:.1f}{'':<16}"
        else:
            row += f"{'N/A':<20}"
    print(row)
    
    print("="*90 + "\n")


def auto_configure_resources(max_nodes: Optional[int] = None, 
                              partition: Optional[str] = None) -> Dict[str, Any]:
    """
    Auto-configure resource settings based on system detection.
    
    Args:
        max_nodes: User-specified max nodes (optional)
        partition: Partition to query (optional)
        
    Returns:
        Dictionary with recommended resource configuration
    """
    detector = SystemDetector()
    resources = detector.detect_all(partition)
    
    recommended_nodes = detector.get_recommended_nodes(max_nodes)
    
    config = {
        'max_nodes': recommended_nodes,
        'procs_per_node': resources.cores_per_node,
        'gpus_per_node': resources.gpus_per_node,
        'memory_per_node_gb': resources.memory_per_node_gb,
        'scheduler': resources.scheduler_type,
        'partition': resources.partition_name,  # Return resolved partition name
        'memory_per_node': f"{resources.memory_per_node_gb:.0f}GB",
        'node_sequence': detector.get_scaling_sequence(recommended_nodes)
    }
    
    logger.info("Auto-configured resources:")
    for key, value in config.items():
        logger.info(f"  {key}: {value}")
    
    return config
