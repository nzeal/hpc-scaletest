# core/types.py
from enum import Enum
from typing import NewType

JobID = NewType('JobID', str)
JobStatus = Enum('JobStatus', 'PENDING RUNNING COMPLETED FAILED')
ScalingType = Enum('ScalingType', 'STRONG WEAK')
BuildBackend = Enum('BuildBackend', 'MAKE CMAKE AUTOTOOLS EASYBUILD SPACK')
SchedulerBackend = Enum('SchedulerBackend', 'LOCAL SLURM')
LauncherBackend = Enum('LauncherBackend', 'SRUN MPIRUN')
ModuleBackend = Enum('ModuleBackend', 'NOMOD TMOD TMOD4 LMOD')
