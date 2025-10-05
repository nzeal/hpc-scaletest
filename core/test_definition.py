# core/test_definition.py
from pathlib import Path
from typing import List, Dict
from .config import TestConfig, ResourceConfig, ScalingConfig
from .types import ScalingType, SchedulerBackend, LauncherBackend, ModuleBackend, BuildBackend
from .factory import BackendFactory

class Test:
    def __init__(self, name: str, input_file: Path, command: List[str]):
        self.config = TestConfig(name=name, input_file=input_file, command=command)

    def set_backend(self, scheduler: str, launcher: str = None, module_system: str = None, build_system: str = None):
        self.config.scheduler = SchedulerBackend[scheduler.upper()]
        if launcher:
            self.config.launcher = LauncherBackend[launcher.upper()]
        if module_system:
            self.config.module_system = ModuleBackend[module_system.upper()]
        if build_system:
            self.config.build_system = BuildBackend[build_system.upper()]

    def set_resources(self, max_nodes: int, procs_per_node: int, gpus_per_node: int = 0, memory_per_node: str = "100GB",
                      time_limit: str = "01:00:00", partition: str = "", account: str = ""):
        self.config.resources = ResourceConfig(
            max_nodes=max_nodes, procs_per_node=procs_per_node, gpus_per_node=gpus_per_node,
            memory_per_node=memory_per_node, time_limit=time_limit, partition=partition, account=account
        )

    def set_scaling(self, scaling_type: str, max_nodes: int, initial_procs: tuple = (1,1,1),
                    initial_domain: tuple = (10.0,10.0,10.0), initial_cells: tuple = (256,256,256)):
        self.config.scaling = ScalingConfig(
            scaling_type=ScalingType[scaling_type.upper()], max_nodes=max_nodes,
            initial_procs=initial_procs, initial_domain=initial_domain, initial_cells=initial_cells
        )

    def set_modules(self, modules: List[str]):
        self.config.modules = modules

    def set_env(self, env: Dict[str, str]):
        self.config.env = env

    def get_config(self) -> TestConfig:
        return self.config
