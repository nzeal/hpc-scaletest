# HPC-ScaleTest Implementation Plan

## Phase 1: Core Interface Alignment (Immediate)

### 1.1 Update BuildSystemInterface

Add missing methods to align interface with concrete implementations:

```python
# core/abstracts.py - ADDITIONS

class BuildSystemInterface(ABC):
    """Abstract interface for build systems."""
    
    # ... existing methods ...
    
    @abstractmethod
    def find_executable(
        self,
        build_dir: Path,
        source_dir: Optional[Path] = None
    ) -> Optional[Path]:
        """
        Find compiled executable in build directory.
        
        Args:
            build_dir: Path to build directory
            source_dir: Optional source directory for name hints
            
        Returns:
            Path to executable or None if not found
        """
        pass
    
    @abstractmethod
    def build_and_find(
        self,
        source_dir: Path,
        build_dir: Path,
        flags: Optional[Dict[str, str]] = None,
        parallel_jobs: int = 4
    ) -> Optional[Path]:
        """
        Complete build workflow with automatic binary detection.
        
        Args:
            source_dir: Path to source directory
            build_dir: Path to build directory
            flags: Build system flags
            parallel_jobs: Number of parallel jobs
            
        Returns:
            Path to compiled executable or None if failed
        """
        pass
```

### 1.2 Add Version Control Interface

```python
# core/abstracts.py - NEW INTERFACE

class VersionControlInterface(ABC):
    """Abstract interface for version control systems."""
    
    def __init__(self, options: Optional[Dict] = None):
        self.options = options or {}
    
    @abstractmethod
    def clone(
        self,
        url: str,
        target_dir: Path,
        branch: Optional[str] = None,
        tag: Optional[str] = None,
        commit: Optional[str] = None,
        shallow: bool = False,
        recursive: bool = True
    ) -> Path:
        """
        Clone repository to target directory.
        
        Args:
            url: Repository URL
            target_dir: Local target directory
            branch: Optional branch to checkout
            tag: Optional tag to checkout
            commit: Optional commit to checkout
            shallow: Whether to do shallow clone
            recursive: Whether to clone submodules
            
        Returns:
            Path to cloned repository
        """
        pass
    
    @abstractmethod
    def checkout(self, repo_dir: Path, ref: str) -> bool:
        """Checkout specific reference (branch/tag/commit)."""
        pass
    
    @abstractmethod
    def get_revision(self, repo_dir: Path) -> str:
        """Get current revision identifier."""
        pass
    
    @abstractmethod
    def pull(self, repo_dir: Path) -> bool:
        """Pull latest changes."""
        pass
```

---

## Phase 2: Plugin Registry Improvements

### 2.1 Enhanced Registry with Build Systems

```python
# core/registry.py - COMPLETE REWRITE

import logging
from typing import Dict, Type, Optional, List, Callable, Any
from dataclasses import dataclass, field
from enum import Enum

logger = logging.getLogger(__name__)


class PluginCategory(Enum):
    """Categories of plugins."""
    SCHEDULER = "scheduler"
    LAUNCHER = "launcher"
    BUILD_SYSTEM = "build_system"
    MODULE_SYSTEM = "module_system"
    VCS = "vcs"
    RESULT_PARSER = "result_parser"
    DIRECTIVE = "directive"
    VALIDATOR = "validator"


@dataclass
class PluginEntry:
    """Registry entry for a plugin."""
    name: str
    cls: Type
    category: PluginCategory
    version: str = "1.0.0"
    description: str = ""
    metadata: Dict = field(default_factory=dict)


class PluginRegistry:
    """
    Central registry for all plugin types.
    
    Thread-safe, instance-based registry with lifecycle management.
    """
    
    _instance: Optional['PluginRegistry'] = None
    _lock = threading.Lock()
    
    def __init__(self):
        """Initialize empty registry."""
        self._plugins: Dict[PluginCategory, Dict[str, PluginEntry]] = {
            category: {} for category in PluginCategory
        }
    
    @classmethod
    def instance(cls) -> 'PluginRegistry':
        """Get singleton instance (thread-safe)."""
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = cls()
        return cls._instance
    
    @classmethod
    def reset(cls):
        """Reset registry (for testing)."""
        with cls._lock:
            cls._instance = None
    
    def register(
        self,
        category: PluginCategory,
        name: str,
        cls: Type,
        version: str = "1.0.0",
        description: str = "",
        metadata: Dict = None
    ):
        """Register a plugin."""
        entry = PluginEntry(
            name=name,
            cls=cls,
            category=category,
            version=version,
            description=description,
            metadata=metadata or {}
        )
        
        if name in self._plugins[category]:
            logger.warning(f"Overwriting {category.value} plugin: {name}")
        
        self._plugins[category][name] = entry
        logger.debug(f"Registered {category.value}: {name} v{version}")
    
    def unregister(self, category: PluginCategory, name: str) -> bool:
        """Unregister a plugin."""
        if name in self._plugins[category]:
            del self._plugins[category][name]
            logger.debug(f"Unregistered {category.value}: {name}")
            return True
        return False
    
    def get(self, category: PluginCategory, name: str) -> Optional[Type]:
        """Get plugin class by category and name."""
        entry = self._plugins[category].get(name)
        return entry.cls if entry else None
    
    def get_entry(self, category: PluginCategory, name: str) -> Optional[PluginEntry]:
        """Get full plugin entry."""
        return self._plugins[category].get(name)
    
    def list(self, category: PluginCategory) -> List[str]:
        """List all plugins in a category."""
        return list(self._plugins[category].keys())
    
    def list_all(self) -> Dict[str, List[str]]:
        """List all registered plugins."""
        return {
            cat.value: list(plugins.keys())
            for cat, plugins in self._plugins.items()
        }
    
    def create_instance(
        self,
        category: PluginCategory,
        name: str,
        options: Dict = None
    ) -> Any:
        """Create plugin instance."""
        cls = self.get(category, name)
        if cls is None:
            available = self.list(category)
            raise ValueError(
                f"Unknown {category.value}: '{name}'. "
                f"Available: {available}"
            )
        return cls(options)


# Decorator for plugin registration
def register_plugin(category: PluginCategory, name: str, **kwargs):
    """Decorator to register a plugin class."""
    def decorator(cls):
        PluginRegistry.instance().register(category, name, cls, **kwargs)
        return cls
    return decorator


# Convenience decorators
def register_scheduler(name: str, **kwargs):
    return register_plugin(PluginCategory.SCHEDULER, name, **kwargs)

def register_launcher(name: str, **kwargs):
    return register_plugin(PluginCategory.LAUNCHER, name, **kwargs)

def register_build_system(name: str, **kwargs):
    return register_plugin(PluginCategory.BUILD_SYSTEM, name, **kwargs)

def register_vcs(name: str, **kwargs):
    return register_plugin(PluginCategory.VCS, name, **kwargs)
```

### 2.2 Updated Factory

```python
# core/factory.py - SIMPLIFIED

from .registry import PluginRegistry, PluginCategory

class BackendFactory:
    """Factory using unified registry."""
    
    @staticmethod
    def create_scheduler(backend: str, options: Dict = None):
        return PluginRegistry.instance().create_instance(
            PluginCategory.SCHEDULER, backend, options
        )
    
    @staticmethod
    def create_launcher(backend: str, options: Dict = None):
        return PluginRegistry.instance().create_instance(
            PluginCategory.LAUNCHER, backend, options
        )
    
    @staticmethod
    def create_build_system(backend: str, options: Dict = None):
        return PluginRegistry.instance().create_instance(
            PluginCategory.BUILD_SYSTEM, backend, options
        )
    
    @staticmethod
    def create_module_system(backend: str, options: Dict = None):
        return PluginRegistry.instance().create_instance(
            PluginCategory.MODULE_SYSTEM, backend, options
        )
    
    @staticmethod
    def create_vcs(backend: str, options: Dict = None):
        return PluginRegistry.instance().create_instance(
            PluginCategory.VCS, backend, options
        )
```

---

## Phase 3: Generic Scaling Engine

### 3.1 Application-Agnostic Scaling

```python
# engine/scaling/generic.py - NEW FILE

"""
Generic scaling engine that uses user-defined variable mappings.
No hardcoded application-specific parameter names.
"""

import copy
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class VariableMapping:
    """
    Mapping of scaling concepts to application-specific variable names.
    
    Example for iPIC3D:
        domain_size: {x: "Lx", y: "Ly", z: "Lz"}
        grid_cells: {x: "nxc", y: "nyc", z: "nzc"}
        mpi_decomposition: {x: "XLEN", y: "YLEN", z: "ZLEN"}
    
    Example for OpenFOAM:
        domain_size: {x: "xMax", y: "yMax", z: "zMax"}
        grid_cells: {x: "nx", y: "ny", z: "nz"}
        mpi_decomposition: {x: "numberOfSubdomains_x", ...}
    """
    domain_size: Dict[str, str] = field(default_factory=dict)
    grid_cells: Dict[str, str] = field(default_factory=dict)
    mpi_decomposition: Dict[str, str] = field(default_factory=dict)
    particles_per_cell: Dict[str, str] = field(default_factory=dict)
    custom: Dict[str, str] = field(default_factory=dict)


@dataclass
class ScalingConfiguration:
    """Configuration for a single scaling point."""
    nodes: int
    total_ranks: int
    mpi_decomposition: Tuple[int, int, int]
    grid_cells: Tuple[int, int, int]
    domain_size: Optional[Tuple[float, float, float]] = None
    custom_values: Dict[str, Any] = field(default_factory=dict)


class GenericScalingEngine:
    """
    Application-agnostic scaling engine.
    
    Uses variable mappings to translate between generic scaling
    concepts and application-specific parameter names.
    """
    
    def __init__(
        self,
        variable_mapping: VariableMapping,
        scaling_type: str = "weak",
        dimensions: int = 3
    ):
        """
        Initialize scaling engine.
        
        Args:
            variable_mapping: Mapping of concepts to variable names
            scaling_type: "weak" or "strong"
            dimensions: 1, 2, or 3 dimensional scaling
        """
        self.var_map = variable_mapping
        self.scaling_type = scaling_type.lower()
        self.dims = dimensions
        
        if self.scaling_type not in ("weak", "strong"):
            raise ValueError(f"Invalid scaling type: {scaling_type}")
        
        if self.dims not in (1, 2, 3):
            raise ValueError(f"Invalid dimensions: {dimensions}")
    
    def generate_configurations(
        self,
        node_sequence: List[int],
        procs_per_node: int,
        base_config: Dict[str, Any]
    ) -> List[ScalingConfiguration]:
        """
        Generate scaling configurations for all nodes in sequence.
        
        Args:
            node_sequence: List of node counts [1, 2, 4, 8, ...]
            procs_per_node: MPI ranks per node
            base_config: Initial configuration values
            
        Returns:
            List of ScalingConfiguration objects
        """
        configs = []
        
        for idx, num_nodes in enumerate(node_sequence):
            total_ranks = num_nodes * procs_per_node
            
            if idx == 0:
                # Baseline configuration
                config = self._create_baseline(base_config, total_ranks)
            elif self.scaling_type == "weak":
                config = self._weak_scale(configs[-1], num_nodes, procs_per_node)
            else:
                config = self._strong_scale(configs[0], num_nodes, procs_per_node)
            
            config.nodes = num_nodes
            config.total_ranks = total_ranks
            configs.append(config)
            
            logger.info(f"Node {num_nodes}: {config}")
        
        return configs
    
    def _create_baseline(
        self,
        base_config: Dict,
        total_ranks: int
    ) -> ScalingConfiguration:
        """Create baseline configuration from user input."""
        # Extract values using variable mapping
        grid = (
            base_config.get('grid_cells', {}).get('x', 64),
            base_config.get('grid_cells', {}).get('y', 64),
            base_config.get('grid_cells', {}).get('z', 64)
        )
        
        domain = None
        if 'domain_size' in base_config:
            domain = (
                base_config['domain_size'].get('x', 1.0),
                base_config['domain_size'].get('y', 1.0),
                base_config['domain_size'].get('z', 1.0)
            )
        
        # Calculate initial MPI decomposition
        mpi_decomp = self._factorize_ranks(total_ranks, self.dims)
        
        return ScalingConfiguration(
            nodes=1,
            total_ranks=total_ranks,
            mpi_decomposition=mpi_decomp,
            grid_cells=grid,
            domain_size=domain,
            custom_values=base_config.get('custom', {})
        )
    
    def _weak_scale(
        self,
        prev_config: ScalingConfiguration,
        num_nodes: int,
        procs_per_node: int
    ) -> ScalingConfiguration:
        """Apply weak scaling - increase problem size with resources."""
        total_ranks = num_nodes * procs_per_node
        
        # Calculate scale factor per dimension
        scale_per_dim = self._calculate_scale_factor(
            prev_config.total_ranks,
            total_ranks
        )
        
        # Scale grid
        new_grid = self._scale_tuple(prev_config.grid_cells, scale_per_dim)
        
        # Scale domain
        new_domain = None
        if prev_config.domain_size:
            new_domain = self._scale_tuple(prev_config.domain_size, scale_per_dim)
        
        # Recalculate MPI decomposition
        new_decomp = self._factorize_ranks(total_ranks, self.dims)
        
        return ScalingConfiguration(
            nodes=num_nodes,
            total_ranks=total_ranks,
            mpi_decomposition=new_decomp,
            grid_cells=new_grid,
            domain_size=new_domain,
            custom_values=prev_config.custom_values.copy()
        )
    
    def _strong_scale(
        self,
        baseline: ScalingConfiguration,
        num_nodes: int,
        procs_per_node: int
    ) -> ScalingConfiguration:
        """Apply strong scaling - fixed problem size, more resources."""
        total_ranks = num_nodes * procs_per_node
        
        # Grid and domain stay the same
        new_decomp = self._factorize_ranks(total_ranks, self.dims)
        
        return ScalingConfiguration(
            nodes=num_nodes,
            total_ranks=total_ranks,
            mpi_decomposition=new_decomp,
            grid_cells=baseline.grid_cells,  # Unchanged
            domain_size=baseline.domain_size,  # Unchanged
            custom_values=baseline.custom_values.copy()
        )
    
    def _factorize_ranks(self, ranks: int, dims: int) -> Tuple[int, ...]:
        """Factorize rank count into dimensions."""
        if dims == 1:
            return (ranks, 1, 1)
        elif dims == 2:
            # Find factors closest to square
            sqrt = int(ranks ** 0.5)
            for i in range(sqrt, 0, -1):
                if ranks % i == 0:
                    return (i, ranks // i, 1)
            return (ranks, 1, 1)
        else:
            # 3D factorization
            cbrt = int(ranks ** (1/3))
            for i in range(cbrt, 0, -1):
                if ranks % i == 0:
                    remaining = ranks // i
                    sqrt = int(remaining ** 0.5)
                    for j in range(sqrt, 0, -1):
                        if remaining % j == 0:
                            return (i, j, remaining // j)
            return (ranks, 1, 1)
    
    def _calculate_scale_factor(
        self,
        prev_ranks: int,
        new_ranks: int
    ) -> Tuple[float, ...]:
        """Calculate scale factor per dimension."""
        ratio = new_ranks / prev_ranks
        
        if self.dims == 1:
            return (ratio, 1.0, 1.0)
        elif self.dims == 2:
            sqrt = ratio ** 0.5
            return (sqrt, sqrt, 1.0)
        else:
            cbrt = ratio ** (1/3)
            return (cbrt, cbrt, cbrt)
    
    def _scale_tuple(
        self,
        values: Tuple,
        factors: Tuple[float, ...]
    ) -> Tuple:
        """Scale tuple values by factors."""
        scaled = []
        for i, (v, f) in enumerate(zip(values, factors)):
            if isinstance(v, int):
                scaled.append(int(round(v * f)))
            else:
                scaled.append(v * f)
        return tuple(scaled)
    
    def apply_to_input_file(
        self,
        config: ScalingConfiguration,
        input_parser,
        input_file: Path,
        output_file: Path
    ) -> bool:
        """
        Apply configuration to input file using variable mapping.
        
        Args:
            config: Scaling configuration to apply
            input_parser: Parser for the input file format
            input_file: Source input file
            output_file: Destination for modified file
            
        Returns:
            True if successful
        """
        # Build modifications dictionary using variable mapping
        modifications = {}
        
        # Grid cells
        if self.var_map.grid_cells:
            for dim, var_name in self.var_map.grid_cells.items():
                idx = {'x': 0, 'y': 1, 'z': 2}.get(dim, 0)
                if idx < len(config.grid_cells):
                    modifications[var_name] = config.grid_cells[idx]
        
        # Domain size
        if self.var_map.domain_size and config.domain_size:
            for dim, var_name in self.var_map.domain_size.items():
                idx = {'x': 0, 'y': 1, 'z': 2}.get(dim, 0)
                if idx < len(config.domain_size):
                    modifications[var_name] = config.domain_size[idx]
        
        # MPI decomposition
        if self.var_map.mpi_decomposition:
            for dim, var_name in self.var_map.mpi_decomposition.items():
                idx = {'x': 0, 'y': 1, 'z': 2}.get(dim, 0)
                if idx < len(config.mpi_decomposition):
                    modifications[var_name] = config.mpi_decomposition[idx]
        
        # Custom values
        for var_name, value in config.custom_values.items():
            if var_name in self.var_map.custom:
                modifications[self.var_map.custom[var_name]] = value
        
        logger.info(f"Applying modifications: {modifications}")
        
        return input_parser.modify_input_file(
            input_file, output_file, modifications
        )
```

---

## Phase 4: Enhanced Code Acquisition

### 4.1 Git Backend Implementation

```python
# plugins/vcs/git.py - NEW FILE

import subprocess
import logging
from pathlib import Path
from typing import Optional

from core.abstracts import VersionControlInterface
from core.registry import register_vcs

logger = logging.getLogger(__name__)


@register_vcs('git', description="Git version control")
class GitBackend(VersionControlInterface):
    """Git version control backend with full feature support."""
    
    def clone(
        self,
        url: str,
        target_dir: Path,
        branch: Optional[str] = None,
        tag: Optional[str] = None,
        commit: Optional[str] = None,
        shallow: bool = False,
        recursive: bool = True
    ) -> Path:
        """Clone repository with options."""
        cmd = ["git", "clone"]
        
        # Shallow clone
        if shallow:
            cmd.extend(["--depth", "1"])
        
        # Branch or tag
        if branch:
            cmd.extend(["--branch", branch])
        elif tag:
            cmd.extend(["--branch", tag])
        
        # Submodules
        if recursive:
            cmd.append("--recursive")
        
        cmd.extend([url, str(target_dir)])
        
        logger.info(f"Cloning: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"Git clone failed: {result.stderr}")
        
        # Checkout specific commit if requested
        if commit:
            self.checkout(target_dir, commit)
        
        logger.info(f"Cloned to {target_dir}")
        return target_dir
    
    def checkout(self, repo_dir: Path, ref: str) -> bool:
        """Checkout specific reference."""
        result = subprocess.run(
            ["git", "checkout", ref],
            cwd=repo_dir,
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            logger.error(f"Checkout failed: {result.stderr}")
            return False
        
        logger.info(f"Checked out: {ref}")
        return True
    
    def get_revision(self, repo_dir: Path) -> str:
        """Get current commit SHA."""
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_dir,
            capture_output=True,
            text=True
        )
        return result.stdout.strip() if result.returncode == 0 else ""
    
    def pull(self, repo_dir: Path) -> bool:
        """Pull latest changes."""
        result = subprocess.run(
            ["git", "pull"],
            cwd=repo_dir,
            capture_output=True,
            text=True
        )
        return result.returncode == 0
    
    def get_remote_url(self, repo_dir: Path) -> str:
        """Get remote origin URL."""
        result = subprocess.run(
            ["git", "remote", "get-url", "origin"],
            cwd=repo_dir,
            capture_output=True,
            text=True
        )
        return result.stdout.strip() if result.returncode == 0 else ""
```

---

## Phase 5: PBS Scheduler Implementation

```python
# plugins/schedulers/pbs.py - NEW FILE

import subprocess
import logging
import re
import time
from pathlib import Path
from typing import List, Optional

from core.abstracts import SchedulerInterface
from core.config import JobConfig, ResourceConfig
from core.types import JobStatus
from core.registry import register_scheduler

logger = logging.getLogger(__name__)


@register_scheduler('pbs', description="PBS/Torque scheduler")
class PBSScheduler(SchedulerInterface):
    """PBS/Torque workload manager backend."""
    
    def generate_job_script(
        self,
        job_config: JobConfig,
        resource_config: ResourceConfig,
        command: List[str],
        env_setup: List[str]
    ) -> str:
        """Generate PBS batch script."""
        
        script = "#!/bin/bash\n\n"
        
        # Resource specification
        ppn = resource_config.procs_per_node
        script += f"#PBS -l nodes={job_config.num_nodes}:ppn={ppn}\n"
        script += f"#PBS -l walltime={resource_config.time_limit}\n"
        
        # Queue/partition
        if resource_config.partition:
            script += f"#PBS -q {resource_config.partition}\n"
        
        # Account
        if resource_config.account:
            script += f"#PBS -A {resource_config.account}\n"
        
        # Job name
        script += f"#PBS -N {job_config.job_id}\n"
        
        # Output files
        script += f"#PBS -o {job_config.working_dir}/job.out\n"
        script += f"#PBS -e {job_config.working_dir}/job.err\n"
        
        # Email notifications
        if resource_config.mail_user:
            script += "#PBS -m abe\n"
            script += f"#PBS -M {resource_config.mail_user}\n"
        
        # GPU resources
        if resource_config.gpus_per_node > 0:
            script += f"#PBS -l gpus={resource_config.gpus_per_node}\n"
        
        # Memory
        if resource_config.memory_per_node:
            script += f"#PBS -l mem={resource_config.memory_per_node}\n"
        
        script += "\n# Environment setup\n"
        for cmd in env_setup:
            script += f"{cmd}\n"
        
        script += "\n# Change to working directory\n"
        script += "cd $PBS_O_WORKDIR\n"
        
        script += "\n# Set OpenMP threads\n"
        script += "export OMP_NUM_THREADS=${PBS_NUM_PPN:-1}\n"
        
        script += "\n# Run application\n"
        total_procs = job_config.num_nodes * ppn
        script += f"mpirun -np {total_procs} {' '.join(command)}\n"
        
        return script
    
    def submit_job(self, script_path: Path) -> str:
        """Submit job and return job ID."""
        result = subprocess.run(
            ["qsub", str(script_path)],
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            raise RuntimeError(f"qsub failed: {result.stderr}")
        
        # Parse job ID from output (format varies by PBS version)
        # Typical format: "12345.hostname" or just "12345"
        job_id = result.stdout.strip().split('.')[0]
        logger.info(f"Submitted job: {job_id}")
        return job_id
    
    def get_job_status(self, job_id: str) -> JobStatus:
        """Query job status."""
        result = subprocess.run(
            ["qstat", "-f", job_id],
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            # Job not found - likely completed
            return JobStatus.COMPLETED
        
        # Parse job state
        state_match = re.search(r'job_state\s*=\s*(\w+)', result.stdout)
        if state_match:
            state = state_match.group(1)
            state_map = {
                'Q': JobStatus.PENDING,
                'R': JobStatus.RUNNING,
                'E': JobStatus.RUNNING,  # Exiting
                'H': JobStatus.PENDING,  # Held
                'C': JobStatus.COMPLETED,
                'F': JobStatus.FAILED,
            }
            return state_map.get(state, JobStatus.UNKNOWN)
        
        return JobStatus.UNKNOWN
    
    def cancel_job(self, job_id: str) -> bool:
        """Cancel a job."""
        result = subprocess.run(
            ["qdel", job_id],
            capture_output=True,
            text=True
        )
        return result.returncode == 0
    
    def wait_for_completion(
        self,
        job_id: str,
        timeout: Optional[int] = None
    ) -> JobStatus:
        """Wait for job to complete."""
        start_time = time.time()
        poll_interval = 30  # seconds
        
        while True:
            status = self.get_job_status(job_id)
            
            if status in (JobStatus.COMPLETED, JobStatus.FAILED):
                return status
            
            if timeout and (time.time() - start_time) > timeout:
                logger.warning(f"Timeout waiting for job {job_id}")
                return JobStatus.TIMEOUT
            
            time.sleep(poll_interval)
```

---

## Summary

This implementation plan provides:

1. **Phase 1**: Interface alignment - ensures abstract interfaces match implementations
2. **Phase 2**: Registry improvements - unified plugin management
3. **Phase 3**: Generic scaling - removes application-specific hardcoding
4. **Phase 4**: VCS support - proper version control abstraction
5. **Phase 5**: PBS scheduler - broader HPC system support

Each phase can be implemented incrementally without breaking existing functionality.
