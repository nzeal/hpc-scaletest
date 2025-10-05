# core/factory.py
from typing import Dict, Any
from .types import SchedulerBackend, LauncherBackend, ModuleBackend, BuildBackend
from ..backends.schedulers.local import LocalScheduler
from ..backends.schedulers.slurm import SlurmScheduler
from ..backends.launchers.srun import SrunLauncher
from ..backends.launchers.mpirun import MpiRunLauncher
from ..backends.modules.nomod import NoModBackend
from ..backends.modules.tmod import TModBackend
from ..backends.modules.tmod4 import TMod4Backend
from ..backends.modules.lmod import LModBackend
from ..backends.builds.make import MakeBackend
from ..backends.builds.cmake import CMakeBackend
from ..backends.builds.autotools import AutotoolsBackend
from ..backends.builds.easybuild import EasyBuildBackend
from ..backends.builds.spack import SpackBackend
from .abstracts import AbstractScheduler, AbstractLauncher, AbstractModuleSystem, AbstractBuildSystem

class BackendFactory:
    _schedulers: Dict[SchedulerBackend, type[AbstractScheduler]] = {
        SchedulerBackend.LOCAL: LocalScheduler,
        SchedulerBackend.SLURM: SlurmScheduler,
    }

    _launchers: Dict[LauncherBackend, type[AbstractLauncher]] = {
        LauncherBackend.SRUN: SrunLauncher,
        LauncherBackend.MPIRUN: MpiRunLauncher,
    }

    _modules: Dict[ModuleBackend, type[AbstractModuleSystem]] = {
        ModuleBackend.NOMOD: NoModBackend,
        ModuleBackend.TMOD: TModBackend,
        ModuleBackend.TMOD4: TMod4Backend,
        ModuleBackend.LMOD: LModBackend,
    }

    _builds: Dict[BuildBackend, type[AbstractBuildSystem]] = {
        BuildBackend.MAKE: MakeBackend,
        BuildBackend.CMAKE: CMakeBackend,
        BuildBackend.AUTOTOOLS: AutotoolsBackend,
        BuildBackend.EASYBUILD: EasyBuildBackend,
        BuildBackend.SPACK: SpackBackend,
    }

    @classmethod
    def create_scheduler(cls, backend: SchedulerBackend, **kwargs) -> AbstractScheduler:
        return cls._schedulers[backend](**kwargs)

    @classmethod
    def create_launcher(cls, backend: LauncherBackend, **kwargs) -> AbstractLauncher:
        return cls._launchers[backend](**kwargs)

    @classmethod
    def create_module_system(cls, backend: ModuleBackend, **kwargs) -> AbstractModuleSystem:
        return cls._modules[backend](**kwargs)

    @classmethod
    def create_build_system(cls, backend: BuildBackend, **kwargs) -> AbstractBuildSystem:
        return cls._builds[backend](**kwargs)
