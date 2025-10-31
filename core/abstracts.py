"""
Abstract base classes defining the plugin interfaces.
"""

import logging
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Dict, Optional, Tuple

from .config import JobConfig, ResourceConfig
from .types import JobStatus


logger = logging.getLogger(__name__)


class SchedulerInterface(ABC):
    """Abstract interface for job schedulers."""
    
    def __init__(self, options: Optional[Dict] = None):
        self.options = options or {}
    
    @abstractmethod
    def generate_job_script(
        self,
        job_config: JobConfig,
        resource_config: ResourceConfig,
        command: List[str],
        env_setup: List[str]
    ) -> str:
        """Generate a job submission script."""
        pass
    
    @abstractmethod
    def submit_job(self, script_path: Path) -> str:
        """Submit a job and return job ID."""
        pass
    
    @abstractmethod
    def get_job_status(self, job_id: str) -> JobStatus:
        """Query the status of a job."""
        pass
    
    @abstractmethod
    def cancel_job(self, job_id: str) -> bool:
        """Cancel a running job."""
        pass
    
    @abstractmethod
    def wait_for_completion(
        self,
        job_id: str,
        timeout: Optional[int] = None
    ) -> JobStatus:
        """Wait for job completion."""
        pass


class LauncherInterface(ABC):
    """Abstract interface for MPI launchers."""
    
    def __init__(self, options: Optional[Dict] = None):
        self.options = options or {}
    
    @abstractmethod
    def generate_launch_command(
        self,
        job_config: JobConfig,
        executable: List[str],
        resource_config: ResourceConfig
    ) -> List[str]:
        """Generate the MPI launch command."""
        pass
    
    @abstractmethod
    def supports_gpu_binding(self) -> bool:
        """Check if launcher supports GPU binding."""
        pass


class ModuleSystemInterface(ABC):
    """Abstract interface for environment module systems."""
    
    def __init__(self, options: Optional[Dict] = None):
        self.options = options or {}
    
    @abstractmethod
    def generate_load_commands(self, modules: List[str]) -> List[str]:
        """Generate commands to load modules."""
        pass
    
    @abstractmethod
    def generate_unload_commands(self, modules: List[str]) -> List[str]:
        """Generate commands to unload modules."""
        pass
    
    @abstractmethod
    def list_available_modules(self, pattern: Optional[str] = None) -> List[str]:
        """List available modules."""
        pass
    
    @abstractmethod
    def is_module_available(self, module: str) -> bool:
        """Check if a module is available."""
        pass


class BuildSystemInterface(ABC):
    """Abstract interface for build systems."""
    
    def __init__(self, options: Optional[Dict] = None):
        self.options = options or {}
    
    @abstractmethod
    def configure(
        self,
        source_dir: Path,
        build_dir: Path,
        flags: Optional[Dict[str, str]] = None
    ) -> bool:
        """Configure the build system."""
        pass
    
    @abstractmethod
    def build(
        self,
        build_dir: Path,
        parallel_jobs: int = 1
    ) -> bool:
        """Execute the build."""
        pass
    
    @abstractmethod
    def install(
        self,
        build_dir: Path,
        install_dir: Path
    ) -> bool:
        """Install the built artifacts."""
        pass
    
    @abstractmethod
    def clean(self, build_dir: Path) -> bool:
        """Clean build artifacts."""
        pass


class ResultParserInterface(ABC):
    """Abstract interface for parsing job results."""
    
    @abstractmethod
    def parse_output(self, output_file: Path) -> Dict:
        """Parse job output and extract metrics."""
        pass
    
    @abstractmethod
    def extract_timing(self, output_file: Path) -> Optional[float]:
        """Extract timing information from output."""
        pass
    
    @abstractmethod
    def extract_performance_metrics(self, output_file: Path) -> Dict:
        """Extract performance metrics."""
        pass


class DefaultResultParser(ResultParserInterface):
    """Default result parser for common benchmark outputs."""
    
    def parse_output(self, output_file: Path) -> Dict:
        try:
            content = output_file.read_text()
            metrics = {}
            
            # Extract wall time using multiple patterns
            time_patterns = [
                r'(?:Time|Wall\s+time|Runtime|Elapsed):\s*(\d+(?:\.\d+)?)',
                r'Total time: (\d+(?:\.\d+)?) seconds',
                r'Performance: .* time = (\d+(?:\.\d+)?)',
                r'Total of (\d+(?:\.\d+)?) seconds elapsed for process'
            ]
            
            for pattern in time_patterns:
                match = re.search(pattern, content, re.I | re.M)
                if match:
                    metrics['wall_time'] = float(match.group(1))
                    break
            
            # Add more extractions as needed (e.g., FLOPS)
            flops_match = re.search(r'(?:GFLOPS|Performance):\s*(\d+(?:\.\d+)?)', content, re.I)
            if flops_match:
                metrics['performance_gflops'] = float(flops_match.group(1))
            
            return metrics
        except Exception as e:
            logger.error(f"Failed to parse {output_file}: {e}")
            return {}
    
    def extract_timing(self, output_file: Path) -> Optional[float]:
        metrics = self.parse_output(output_file)
        return metrics.get('wall_time')
    
    def extract_performance_metrics(self, output_file: Path) -> Dict:
        return self.parse_output(output_file)
