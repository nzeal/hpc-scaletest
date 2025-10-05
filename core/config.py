# core/config.py
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
from pathlib import Path
from .types import ScalingType, SchedulerBackend, LauncherBackend, ModuleBackend, BuildBackend

@dataclass
class ResourceConfig:
    max_nodes: int
    procs_per_node: int
    gpus_per_node: int = 0
    memory_per_node: str = "100GB"
    time_limit: str = "01:00:00"
    partition: str = ""
    account: str = ""

@dataclass
class ScalingConfig:
    scaling_type: ScalingType
    max_nodes: int
    initial_procs: Tuple[int, int, int]
    initial_domain: Tuple[float, float, float]
    initial_cells: Tuple[int, int, int]

@dataclass
class TestConfig:
    name: str
    input_file: Path
    command: List[str]
    scheduler: SchedulerBackend = SchedulerBackend.LOCAL
    launcher: LauncherBackend = LauncherBackend.MPIRUN
    module_system: ModuleBackend = ModuleBackend.NOMOD
    build_system: Optional[BuildBackend] = None
    modules: List[str] = None
    env: Dict[str, str] = None
    resources: ResourceConfig = None
    scaling: ScalingConfig = None

    def __post_init__(self):
        if self.modules is None:
            self.modules = []
        if self.env is None:
            self.env = {}
