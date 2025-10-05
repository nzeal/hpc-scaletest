# core/abstracts.py
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional
from pathlib import Path
from .types import JobID, JobStatus

class AbstractScheduler(ABC):
    @abstractmethod
    def submit(self, job_script: Path) -> JobID:
        """Submit a job script and return job ID."""
        pass

    @abstractmethod
    def monitor(self, job_id: JobID) -> JobStatus:
        """Monitor job status."""
        pass

    @abstractmethod
    def cancel(self, job_id: JobID) -> bool:
        """Cancel a job."""
        pass

    @abstractmethod
    def wait(self, job_id: JobID, timeout: Optional[int] = None) -> Dict[str, Any]:
        """Wait for job completion and return results."""
        pass

class AbstractLauncher(ABC):
    @abstractmethod
    def launch(self, command: List[str], nodes: int, procs_per_node: int, env: Dict[str, str] = None) -> Dict[str, Any]:
        """Launch a parallel job and return results (e.g., timings)."""
        pass

class AbstractModuleSystem(ABC):
    @abstractmethod
    def load(self, modules: List[str]) -> None:
        """Load modules."""
        pass

    @abstractmethod
    def unload(self, modules: List[str]) -> None:
        """Unload modules."""
        pass

    @abstractmethod
    def purge(self) -> None:
        """Purge all modules."""
        pass

    @abstractmethod
    def list(self) -> List[str]:
        """List loaded modules."""
        pass

class AbstractBuildSystem(ABC):
    @abstractmethod
    def build(self, source_dir: Path, flags: Dict[str, str] = None) -> Path:
        """Build the project and return build path."""
        pass

    @abstractmethod
    def install(self, build_path: Path) -> Path:
        """Install the build and return install path."""
        pass
