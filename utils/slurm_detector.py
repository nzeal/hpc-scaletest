"""
Automatic SLURM feature detection and partition selection.
Queries Leonardo's Slurm configuration to detect:
  - Available partitions and their limits (MaxNodes, MaxTime, etc.)
  - Per-partition GPU availability
  - QOS options and their specifications
  - Node hardware (cores, memory, CPU model)
  - Per-user limits and scheduling constraints
"""

import subprocess
import logging
import re
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class SlurmPartitionInfo:
    """Information about a SLURM partition."""
    name: str
    total_nodes: int
    available_nodes: int
    max_nodes: int
    max_cpus_per_node: int
    max_time_minutes: int
    allow_accounts: List[str] = field(default_factory=list)
    deny_accounts: List[str] = field(default_factory=list)
    gpu_per_node: int = 0
    has_gres_tmpfs: bool = False
    priority: int = 0  # Higher = better priority for large jobs
    
    def can_run_job(self, num_nodes: int, num_cpus_per_node: int, walltime_minutes: int) -> Tuple[bool, str]:
        """Check if partition can satisfy job requirements."""
        if num_nodes > self.max_nodes:
            return False, f"Requested {num_nodes} nodes exceeds max {self.max_nodes}"
        
        if num_cpus_per_node > self.max_cpus_per_node:
            return False, f"Requested {num_cpus_per_node} CPUs exceeds max {self.max_cpus_per_node}/node"
        
        if walltime_minutes > self.max_time_minutes:
            return False, f"Requested {walltime_minutes}min exceeds max {self.max_time_minutes}min"
        
        # NOTE: Use max_nodes for feasibility check, not available_nodes
        # available_nodes from sinfo only shows currently free nodes, not total partition capacity
        # max_nodes from scontrol is authoritative for what the partition can support
        # We check max_nodes above already, so if we get here, the job is feasible
        
        return True, "OK"
    
    def score_for_job(self, num_nodes: int) -> int:
        """Score partition for job fit (higher is better for large jobs)."""
        # Prefer partitions that can comfortably fit the requested nodes
        # with room to spare for other users
        score = 0
        
        # Primary factor: fit for job size (most important)
        if num_nodes <= self.max_nodes:
            # For large jobs, prefer dedicated high-capacity partitions over unlimited ones
            if num_nodes > 16:
                if 32 <= num_nodes <= 64 and self.max_nodes >= 64:
                    # Perfect fit for 32-64 range (e.g., boost_usr_prod)
                    score += 20000
                elif num_nodes > 16 and self.max_nodes >= 128:
                    score += 10000  # Excellent for large jobs
                elif num_nodes > 16 and self.max_nodes >= 64:
                    score += 9000   # Good for large jobs
                elif num_nodes > 16 and self.max_nodes >= 32:
                    score += 1000   # OK for large jobs
                else:
                    score -= 10000  # Not suitable (too small for this job)
            else:
                # For small jobs (<=16 nodes), prefer smaller specialized partitions
                # This avoids blocking large-job partitions
                if self.max_nodes <= 16:
                    score += 1000   # Perfect fit for small partition
                elif self.max_nodes <= 64:
                    score += 500    # OK, but not ideal (ties up larger partition)
                else:
                    score += 100    # Not ideal (very large partition wasted on small job)
        
        # Secondary factor: availability (very small influence)
        score += self.available_nodes
        
        # Tertiary factor: explicit priority (minimal influence)
        score += self.priority
        
        return score


@dataclass
class SlurmNodeInfo:
    """Information about SLURM compute node."""
    name: str
    cpus: int
    memory_gb: float
    gpu_count: int
    gpu_model: Optional[str] = None


class SlurmDetector:
    """Detect and query SLURM configuration on Leonardo."""
    
    def __init__(self):
        self.partitions: Dict[str, SlurmPartitionInfo] = {}
        self.nodes: Dict[str, SlurmNodeInfo] = {}
        self.qos_list: List[str] = []
        self.user_limits: Dict[str, int] = {}  # user limits like max concurrent jobs
        self._detect_all()
    
    def _run_scontrol(self, args: List[str]) -> str:
        """Run scontrol and return output."""
        try:
            result = subprocess.run(["scontrol"] + args, capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                logger.warning(f"scontrol {' '.join(args)} failed: {result.stderr}")
                return ""
            return result.stdout
        except Exception as e:
            logger.warning(f"Could not run scontrol: {e}")
            return ""
    
    def _run_sinfo(self, args: List[str]) -> str:
        """Run sinfo and return output."""
        try:
            result = subprocess.run(["sinfo"] + args, capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                logger.warning(f"sinfo {' '.join(args)} failed: {result.stderr}")
                return ""
            return result.stdout
        except Exception as e:
            logger.warning(f"Could not run sinfo: {e}")
            return ""
    
    def _detect_all(self):
        """Auto-detect all SLURM features."""
        logger.info("Detecting SLURM configuration...")
        self._detect_partitions()
        self._detect_nodes()
        self._detect_qos()
        self._detect_user_limits()
        logger.info(f"  ✓ Found {len(self.partitions)} partitions")
        logger.info(f"  ✓ Found {len(self.nodes)} node types")
        logger.info(f"  ✓ Found {len(self.qos_list)} QOS options")
    
    def _detect_partitions(self):
        """Detect available partitions and their limits."""
        output = self._run_scontrol(["show", "partition"])
        if not output:
            logger.warning("Could not detect partitions")
            return
        # Parse scontrol output token-by-token (handles multiple key=value per line)
        current_partition: Dict[str, str] = {}
        for raw_line in output.split('\n'):
            line = raw_line.strip()
            # Blank line separates partition blocks in some scontrol formats
            if not line:
                if current_partition:
                    self._parse_partition_info(current_partition)
                    current_partition = {}
                continue

            # Each line can contain multiple Key=Value tokens
            tokens = line.split()
            for tok in tokens:
                if '=' not in tok:
                    continue
                key, value = tok.split('=', 1)
                # If we encounter a new PartitionName and we already have data, flush previous
                if key == 'PartitionName' and current_partition:
                    self._parse_partition_info(current_partition)
                    current_partition = {}
                current_partition[key] = value

        # Flush last partition block if present
        if current_partition:
            self._parse_partition_info(current_partition)
    
    def _parse_partition_info(self, info: Dict[str, str]):
        """Parse partition information from scontrol output."""
        name = info.get('PartitionName', 'unknown')
        
        # Get max nodes (Leonardo sets this per partition)
        max_nodes_str = info.get('MaxNodes', '0')
        # MaxNodes can be "N" or "UNLIMITED"
        if max_nodes_str.upper() == 'UNLIMITED':
            max_nodes = 10000  # Large number for "unlimited"
        else:
            try:
                max_nodes = int(max_nodes_str)
            except ValueError:
                max_nodes = 0
        
        # Get max CPUs per node
        max_cpus = int(info.get('MaxCPUsPerNode', 'UNLIMITED').replace('UNLIMITED', '10000')) or 112
        
        # Parse max time (format: "days-hours:minutes:seconds" or "minutes")
        max_time_str = info.get('MaxTime', '1-00:00:00')
        max_time_minutes = self._parse_time_to_minutes(max_time_str)
        
        # Get GRES info (GPUs, tmpfs, etc.)
        gres_str = info.get('GRES', '')
        gpu_per_node = self._parse_gres_gpus(gres_str)
        has_gres_tmpfs = 'tmpfs' in gres_str.lower()
        
        # Get node count from sinfo
        total_nodes, available_nodes = self._count_partition_nodes(name)
        
        # Priority (lower partition name = higher priority on Leonardo)
        priority = 0
        if 'boost' in name.lower():
            priority = 100
        elif 'dcgp' in name.lower():
            priority = 50
        elif 'gpu' in name.lower():
            priority = 40
        
        partition = SlurmPartitionInfo(
            name=name,
            total_nodes=total_nodes,
            available_nodes=available_nodes,
            max_nodes=max_nodes,
            max_cpus_per_node=max_cpus,
            max_time_minutes=max_time_minutes,
            gpu_per_node=gpu_per_node,
            has_gres_tmpfs=has_gres_tmpfs,
            priority=priority
        )
        
        self.partitions[name] = partition
        logger.debug(f"  Partition {name}: max_nodes={max_nodes}, cpus={max_cpus}, gpus={gpu_per_node}")
    
    def _parse_time_to_minutes(self, time_str: str) -> int:
        """Parse SLURM time format to minutes."""
        if not time_str or time_str.upper() == 'UNLIMITED':
            return 1440 * 365  # 1 year default
        
        total_minutes = 0
        
        # Format: "days-hours:minutes:seconds"
        if '-' in time_str:
            days_part, time_part = time_str.split('-', 1)
            total_minutes += int(days_part) * 1440
            time_str = time_part
        
        # Now parse "hours:minutes:seconds"
        parts = time_str.split(':')
        if len(parts) >= 3:
            total_minutes += int(parts[0]) * 60 + int(parts[1])
        elif len(parts) == 2:
            total_minutes += int(parts[0]) * 60 + int(parts[1])
        elif len(parts) == 1:
            total_minutes += int(parts[0])
        
        return total_minutes
    
    def _parse_gres_gpus(self, gres_str: str) -> int:
        """Extract GPU count from GRES string."""
        # Format: "gpu:type:count" or similar
        match = re.search(r'gpu:(?:\w+:)?(\d+)', gres_str, re.IGNORECASE)
        if match:
            return int(match.group(1))
        return 0
    
    def _count_partition_nodes(self, partition_name: str) -> Tuple[int, int]:
        """Count total and available nodes in partition."""
        output = self._run_sinfo(["-p", partition_name, "-h", "-o", "%t"])
        if not output:
            return 0, 0
        
        total = 0
        available = 0
        
        for line in output.split('\n'):
            state = line.strip()
            if state:
                total += 1
                # Available states: idle, allocated (can't use), mixed, etc.
                if state in ['idle', 'mixed']:
                    available += 1
        
        return total, available
    
    def _detect_nodes(self):
        """Detect node hardware configuration."""
        output = self._run_sinfo(["-N", "-h", "-o", "%N %C %m %G"])
        if not output:
            logger.warning("Could not detect node configuration")
            return
        
        # Parse each line
        for line in output.split('\n'):
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 4:
                name = parts[0]
                cpus_str = parts[1]  # Format: "allocated/idle/other/total"
                memory_str = parts[2]
                gpu_str = parts[3]
                
                # Parse CPUs (take total)
                if cpus_str.count('/') == 3:
                    cpu_total = int(cpus_str.split('/')[-1])
                else:
                    cpu_total = int(cpus_str) if cpus_str.isdigit() else 112  # default
                
                # Parse memory
                memory_gb = self._parse_memory(memory_str)
                
                # Parse GPUs
                gpu_count = 0
                if gpu_str != '(null)':
                    # Format: "gpu:type:count" or "count(type)"
                    match = re.search(r'(\d+)', gpu_str)
                    if match:
                        gpu_count = int(match.group(1))
                
                self.nodes[name] = SlurmNodeInfo(
                    name=name,
                    cpus=cpu_total,
                    memory_gb=memory_gb,
                    gpu_count=gpu_count
                )
    
    def _parse_memory(self, mem_str: str) -> float:
        """Parse memory string to GB."""
        match = re.match(r'(\d+)([KMG])?', mem_str)
        if match:
            value = int(match.group(1))
            unit = match.group(2) or 'M'
            
            if unit == 'K':
                return value / (1024 * 1024)
            elif unit == 'M':
                return value / 1024
            elif unit == 'G':
                return float(value)
        
        return 0.0
    
    def _detect_qos(self):
        """Detect available QOS options."""
        output = self._run_scontrol(["show", "qos"])
        if not output:
            return
        
        for line in output.split('\n'):
            if line.startswith('Name='):
                qos_name = line.split('=', 1)[1]
                self.qos_list.append(qos_name)
    
    def _detect_user_limits(self):
        """Detect per-user limits."""
        try:
            result = subprocess.run(["sshare", "-u", "-n"], capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                # Parse sshare output
                for line in result.stdout.split('\n')[1:]:  # Skip header
                    parts = line.split()
                    if len(parts) >= 2:
                        user = parts[0]
                        # Could extract max job limits here
        except:
            pass
    
    def select_best_partition(self, num_nodes: int, num_cpus_per_node: int,
                             walltime_minutes: int = 120) -> Optional[str]:
        """
        Select the best partition for a job.
        
        Returns partition name that can run the job, or None if impossible.
        """
        candidates = []
        
        for partition_name, partition_info in self.partitions.items():
            can_run, reason = partition_info.can_run_job(num_nodes, num_cpus_per_node, walltime_minutes)
            
            if can_run:
                score = partition_info.score_for_job(num_nodes)
                candidates.append((score, partition_name))
                logger.debug(f"  {partition_name}: score={score}, available={partition_info.available_nodes}/{partition_info.max_nodes}")
            else:
                logger.debug(f"  {partition_name}: cannot run - {reason}")
        
        if not candidates:
            logger.warning(f"No partition can run job: {num_nodes} nodes, {num_cpus_per_node} cpus, {walltime_minutes}min")
            return None
        
        # Sort by score (descending) and return best
        candidates.sort(reverse=True)
        best_partition = candidates[0][1]
        logger.info(f"  ✓ Selected partition: {best_partition}")
        
        return best_partition
    
    def get_node_cores(self, partition: str = None) -> int:
        """Get number of cores per node (auto-detect or from partition)."""
        if not self.nodes:
            return 112  # Leonardo default
        
        # Get first node info
        first_node = next(iter(self.nodes.values()))
        return first_node.cpus
    
    def get_node_memory_gb(self, partition: str = None) -> float:
        """Get memory per node in GB (auto-detect or from partition)."""
        if not self.nodes:
            return 502.0  # Leonardo dcgp default
        
        # Get first node info
        first_node = next(iter(self.nodes.values()))
        return first_node.memory_gb
    
    def get_default_walltime(self) -> str:
        """Get reasonable default walltime for partition."""
        return "02:00:00"  # 2 hours default


# Singleton instance
_detector = None


def get_slurm_detector() -> SlurmDetector:
    """Get or create SLURM detector singleton."""
    global _detector
    if _detector is None:
        _detector = SlurmDetector()
    return _detector
