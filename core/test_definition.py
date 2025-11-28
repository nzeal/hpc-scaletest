"""
User-facing Test definition class.
Simple API: Set name, command (exe + input), source_dir for auto-build, resources, scaling.
"""

from pathlib import Path
from typing import List, Dict, Optional

from .config import (
    BackendConfig, ResourceConfig, ScalingConfig,
    EnvironmentConfig, BuildConfig
)
from .types import (
    ScalingType, SchedulerBackend, LauncherBackend,
    ModuleBackend, BuildBackend,
    ProcsDecomposition, DomainSize, CellCount,
    DEFAULT_OUTPUT_DIR
)


class Test:
    """Main test definition class."""
    
    def __init__(
        self,
        name: str,
        input_file: Optional[Path] = None,
        command: Optional[List[str]] = None,
        output_dir: Path = Path("output")
    ):
        self.name = name
        self.input_file = input_file
        self.command = command or []
        self.output_dir = output_dir
        
        # Configuration objects
        self.backend_config = BackendConfig()
        self.resource_config = ResourceConfig()
        self.scaling_config = ScalingConfig()
        self.environment_config = EnvironmentConfig()
        self.build_config = BuildConfig()
        
        # New attribute to control automatic job submission
        self.auto_submit_jobs = True  # By default, submit jobs automatically
        
        # LLM configuration for intelligent parameter mapping
        self.use_llm_discovery = False
        self.openai_api_key: Optional[str] = None
        self.llm_model = "gpt-4"
        self.parameter_mapping: Optional[Dict[str, List[str]]] = None
    
    def set_backend(
        self,
        scheduler: str = "slurm",
        launcher: str = "mpirun",  # Default to mpirun for simplicity
        module_system: str = "lmod"
    ) -> 'Test':
        """Simple backend setup."""
        self.backend_config.scheduler = SchedulerBackend(scheduler)
        self.backend_config.launcher = LauncherBackend(launcher)
        self.backend_config.module_system = ModuleBackend(module_system)
        self.backend_config.build_system = BuildBackend.CMAKE  # Default
        return self
    
    def set_resources(
        self,
        max_nodes: int = 1,
        procs_per_node: int = 32,
        gpus_per_node: int = 0,
        time_limit: str = "02:00:00",
        partition: str = "X_usr_prod",
        account: str = "cin_X",
        qos: Optional[str] = None,
        qos_mapping: Optional[Dict] = None,
        exclusive: bool = False,
        mail_user: Optional[str] = None
    ) -> 'Test':
        """Simple resource setup."""
        self.resource_config.max_nodes = max_nodes
        self.resource_config.procs_per_node = procs_per_node
        self.resource_config.gpus_per_node = gpus_per_node
        self.resource_config.time_limit = time_limit
        self.resource_config.partition = partition
        self.resource_config.account = account
        self.resource_config.qos = qos
        self.resource_config.qos_mapping = qos_mapping
        self.resource_config.exclusive = exclusive
        self.resource_config.mail_user = mail_user
        return self
    
    def set_scaling(
        self,
        scaling_type: str = "weak",
        max_nodes: int = 4,
        initial_procs: ProcsDecomposition = (1, 1, 1),
        initial_domain: Optional[DomainSize] = None,
        initial_cells: Optional[CellCount] = None
    ) -> 'Test':
        """Simple scaling setup."""
        self.scaling_config.scaling_type = ScalingType(scaling_type)
        self.scaling_config.max_nodes = max_nodes
        self.scaling_config.initial_procs = initial_procs
        self.scaling_config.initial_domain = initial_domain
        self.scaling_config.initial_cells = initial_cells
        return self
    
    def set_modules(self, modules: List[str]) -> 'Test':
        """Set modules (loaded in script)."""
        self.environment_config.modules = modules
        return self
    
    def set_env(self, env_vars: Dict[str, str]) -> 'Test':
        """Set env vars (exported in script)."""
        self.environment_config.env_vars = env_vars
        return self
    
    def get_input_content(self, job_config) -> str:
        """Generate/modify input for job (scale for weak). Override if needed."""
        if self.input_file:
            base_content = self.input_file.read_text()
            lines = base_content.splitlines()
            new_lines = []
            for line in lines:
                stripped = line.strip()
                if stripped.startswith("Lx"):
                    new_lines.append(f"Lx = {job_config.domain_size[0] if job_config.domain_size else '40.96'}")
                elif stripped.startswith("Ly"):
                    new_lines.append(f"Ly = {job_config.domain_size[1] if job_config.domain_size else '20.48'}")
                elif stripped.startswith("Lz"):
                    new_lines.append(f"Lz = {job_config.domain_size[2] if job_config.domain_size else '1.0'}")
                elif stripped.startswith("nxc"):
                    new_lines.append(f"nxc = {int(job_config.cell_count[0]) if job_config.cell_count else 896}")
                elif stripped.startswith("nyc"):
                    new_lines.append(f"nyc = {int(job_config.cell_count[1]) if job_config.cell_count else 512}")
                elif stripped.startswith("nzc"):
                    new_lines.append(f"nzc = {int(job_config.cell_count[2]) if job_config.cell_count else 1}")
                elif stripped.startswith("XLEN"):
                    new_lines.append(f"XLEN = {job_config.procs_decomposition[0]}")
                elif stripped.startswith("YLEN"):
                    new_lines.append(f"YLEN = {job_config.procs_decomposition[1]}")
                elif stripped.startswith("ZLEN"):
                    new_lines.append(f"ZLEN = {job_config.procs_decomposition[2]}")
                else:
                    new_lines.append(line)
            return "\n".join(new_lines) + "\n"
        return ""
    
    def set_auto_submit(self, auto_submit: bool = True) -> 'Test':
        """Control whether jobs are automatically submitted after generation.
        
        Args:
            auto_submit: If True, jobs are submitted automatically. If False, 
                        only job scripts are generated.
        
        Returns:
            self: For method chaining
        """
        self.auto_submit_jobs = auto_submit
        return self
    
    def validate(self) -> bool:
        """Validate config."""
        if not self.name:
            raise ValueError("Test name is required")
        if not self.command:
            raise ValueError("Test command is required")
        if self.scaling_config.max_nodes < 1:
            raise ValueError("max_nodes must be >= 1")
        if self.build_config.source_dir and not self.build_config.source_dir.exists():
            raise ValueError(f"Source dir not found: {self.build_config.source_dir}")
        if self.input_file and not self.input_file.exists():
            raise ValueError(f"Input file not found: {self.input_file}")
        return True
