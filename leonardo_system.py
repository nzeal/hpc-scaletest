"""
LEONARDO HPC System Configuration

This file demonstrates how to define custom launchers and system configurations
for the LEONARDO supercomputer.

Users can create similar files for their own HPC systems.
"""

import json
import re
from typing import List
from core.registry import register_launcher, JobLauncher
from core.config import JobConfig, ResourceConfig


# ============================================================================
# Custom Launcher Definitions
# ============================================================================

@register_launcher('mpirun-mapby')
class MpirunMapbyLauncher(JobLauncher):
    """
    MPI launcher with custom socket mapping and core binding.
    
    This launcher uses mpirun with specific mapping options for
    optimal process placement on NUMA nodes.
    """
    
    def command(self, job: JobConfig, resource_config: ResourceConfig) -> List[str]:
        """Generate mpirun command with custom mapping."""
        num_cpus_per_task = self.options.get('num_cpus_per_task', 
                                              resource_config.procs_per_node // job.num_procs 
                                              if job.num_procs <= resource_config.procs_per_node 
                                              else 1)
        
        return [
            'mpirun', 
            '-np', str(job.num_procs),
            f'--map-by', f'socket:PE={num_cpus_per_task}',
            '--rank-by', 'core',
            '--report-bindings'
        ]


@register_launcher('mpirun-nsys')
class MpirunNsysLauncher(JobLauncher):
    """
    MPI launcher with NVIDIA Nsight Systems profiling.
    
    This launcher wraps the application with nsys for performance profiling,
    generating per-rank trace files.
    """
    
    def command(self, job: JobConfig, resource_config: ResourceConfig) -> List[str]:
        """Generate mpirun command with Nsight profiling."""
        output_prefix = self.options.get('nsys_output_prefix', 'report_rank')
        trace_types = self.options.get('nsys_trace', 'mpi,cuda')
        
        return [
            'mpirun', 
            '-np', str(job.num_procs), 
            '--report-bindings',
            'nsys', 'profile', 
            '-t', trace_types,
            '--force-overwrite=true',
            '-o', f'{output_prefix}%q{{OMPI_COMM_WORLD_RANK}}', 
            '--stats=true'
        ]


@register_launcher('srun-gpu')
class SrunGpuLauncher(JobLauncher):
    """
    Slurm srun launcher with GPU binding support.
    
    Properly binds MPI ranks to GPUs on multi-GPU nodes.
    """
    
    def command(self, job: JobConfig, resource_config: ResourceConfig) -> List[str]:
        """Generate srun command with GPU binding."""
        cmd = [
            'srun',
            '-n', str(job.num_procs),
            '--ntasks-per-node', str(resource_config.procs_per_node)
        ]
        
        if resource_config.gpus_per_node > 0:
            cmd.extend([
                '--gpus-per-node', str(resource_config.gpus_per_node),
                '--gpu-bind', 'closest'
            ])
        
        return cmd
    
    def supports_gpu_binding(self) -> bool:
        return True


# ============================================================================
# System Configuration
# ============================================================================

site_configuration = {
    "systems": [
        {
            "name": "leonardo",
            "descr": "LEONARDO HPC System at CINECA",
            "hostnames": [r"^login\d+\.leonardo\.local$"],
            "modules_system": "tmod4",
            "partitions": [
                {
                    "name": "login",
                    "descr": "LEONARDO Login Nodes",
                    "scheduler": "local",
                    "launcher": "local",
                    "modules": [],
                    "access": [],
                    "max_jobs": 1,
                    "environs": [
                        "default", 
                        "gcc", 
                        "openmpi", 
                    ],
                    "processor": {
                        "arch": "icelake",
                        "platform": "x86_64",
                        "num_cpus": 32,
                        "num_cpus_per_core": 1,
                        "num_cpus_per_socket": 32,
                        "num_sockets": 1,
                    },
                },
                {
                    "name": "booster",
                    "descr": "LEONARDO Booster partition (GPU nodes)",
                    "scheduler": "slurm",
                    "sched_options": {
                        "use_nodes_option": True
                    },
                    "launcher": "srun",
                    "modules": [],
                    "container_platforms": [
                        {
                            "type": "Singularity",
                            "default": True,
                            "modules": [],
                            "env_vars": [['ENV_VAR', 'VALUE']]
                        }
                    ],
                    "access": [
                        "--partition=boost_usr_prod",
                        "--account=cin_staff"  # Default account - users can override with --account flag
                    ],
                    "resources": [
                        {"name": "qos", "options": ["--qos={qos}"]},
                        {"name": "account", "options": ["--account={account}"]},
                        {"name": "gpu", "options": ["--gres=gpu:{num_gpus_per_node}"]},
                        {"name": "cpufreq", "options": ["--cpu-freq={cpufreq}"]},
                    ],
                    "max_jobs": 10,
                    "environs": [
                        "default",
                        "gcc",
                        "nvhpc",
                        "openmpi",
                        "cuda",
                        "openmpi-gcc-cuda",
                    ],
                    "processor": {
                        "arch": "icelake",
                        "platform": "x86_64",
                        "num_cpus": 32,
                        "num_cpus_per_core": 1,
                        "num_cpus_per_socket": 32,
                        "num_sockets": 1,
                    },
                    "extras": {
                        "l2_cache_per_core_kb": 512,
                        "l3_cache_per_socket_kb": 4096
                    },
                    "devices": [
                        {
                            "type": "gpu",
                            "arch": "sm_80",
                            "model": "A100",
                            "num_devices": 4,
                        }
                    ],
                },
                {
                    "name": "dcgp",
                    "descr": "LEONARDO Data Centric General Purpose (DCGP) partition",
                    "scheduler": "slurm",
                    "launcher": "srun",
                    "max_jobs": 8,
                    "modules": [],
                    "access": [
                        "--partition=dcgp_usr_prod",
                        "--account=cin_staff"  # Default account - users can override with --account flag
                    ],
                    "resources": [
                        {"name": "qos", "options": ["--qos={qos}"]},
                        {"name": "account", "options": ["--account={account}"]},
                    ],
                    "environs": ["default", "gcc"],
                    "processor": {
                        "num_cpus": 112,
                        "num_cpus_per_core": 1,
                        "num_cpus_per_socket": 56,
                        "num_sockets": 2,
                    },
                    "extras": {
                        "l2_cache_per_core_kb": 2048,
                        "l3_cache_per_socket_kb": 107520
                    },
                },
            ],
        },
    ],
    "environments": [
        {
            "name": "default",
            "modules": [],
            "cc": "gcc",
            "cxx": "g++",
            "ftn": "gfortran",
            "features": ["openmp"],
            "extras": {"omp_flag": "-fopenmp"},
        },
        {
            "name": "gcc",
            "modules": ["gcc/12.2.0"],
            "cc": "gcc",
            "cxx": "g++",
            "ftn": "gfortran",
            "features": ["openmp"],
            "extras": {"omp_flag": "-fopenmp"},
        },
        {
            "name": "nvhpc",
            "modules": ["nvhpc/24.5"],
            "cc": "nvc",
            "cxx": "nvc++",
            "ftn": "nvfortran",
            "features": ["openmp"],
            "extras": {"omp_flag": "-mp=multicore"},
        },
        {
            "name": "openmpi",
            "modules": ["openmpi/4.1.6--gcc--12.2.0"],
            "cc": "gcc",
            "cxx": "g++",
            "ftn": "gfortran",
            "features": ["openmp", "mpi"],
            "extras": {"omp_flag": "-fopenmp"},
        },
        {
            "name": "openmpi-nvhpc",
            "modules": ["nvhpc/24.5", "openmpi/4.1.6--nvhpc--24.5"],
            "cc": "nvc",
            "cxx": "nvc++",
            "ftn": "nvfortran",
            "features": ["openmp", "mpi", "cuda"],
            "extras": {"omp_flag": "-fopenmp"},
        },
        {
            "name": "cuda",
            "modules": ["cuda/12.3"],
            "cc": "gcc",
            "cxx": "g++",
            "ftn": "gfortran",
            "features": ["openmp", "cuda"],
            "extras": {"omp_flag": "-fopenmp"},
        },
        {
            "name": "openmpi-gcc-cuda",
            "modules": ["openmpi/4.1.6--gcc--12.2.0-cuda-12.2"],
            "cc": "gcc",
            "cxx": "g++",
            "ftn": "gfortran",
            "nvcc": "nvcc",
            "features": ["openmp", "mpi", "cuda"],
            "extras": {"omp_flag": "-fopenmp"},
        },
    ],
    "logging": [
        {
            "level": "debug",
            "handlers": [
                {
                    "type": "file",
                    "name": "reframe.log",
                    "level": "debug",
                    "format": "[%(asctime)s] %(levelname)s: %(check_name)s: %(message)s",
                    "append": False,
                },
                {
                    "type": "stream",
                    "name": "stdout",
                    "level": "info",
                    "format": "%(message)s",
                },
                {
                    "type": "file",
                    "name": "reframe.out",
                    "level": "info",
                    "format": "%(message)s",
                    "append": False,
                },
            ],
            "handlers_perflog": [
                {
                    "type": "filelog",
                    "prefix": "%(check_system)s/%(check_partition)s",
                    "level": "info",
                    "format": (
                        "%(asctime)s|"
                        "reframe %(version)s|"
                        "%(check_job_completion_time)s|"
                        "%(check_info)s|"
                        "%(check_modules)s|"
                        "%(check_result)s|"
                        "%(check_executable)s|"
                        "%(check_executable_opts)s|"
                        "%(check_system)s|"
                        "%(check_partition)s|"
                        "%(check_environ)s|"
                        "%(check_descr)s|"
                        "%(check_job_nodelist)s|"
                        "%(check_num_tasks_per_node)s|"
                        "%(check_num_cpus_per_task)s|"
                        "%(check_num_gpus_per_node)s|"
                        "%(check_num_tasks)s|"
                        "%(check_exclusive_access)s|"
                        "%(check_perfvalues)s"
                    ),
                    "format_perfvars": ("%(check_perf_value)s|" "%(check_perf_unit)s|"),
                    "append": True,
                },
            ],
        }
    ],
    "storage": [
        {
            "enable": True,
            "backend": "postgresql",
            "postgresql_driver": "psycopg2",
            "postgresql_host": "131.175.204.237",
            "postgresql_port": 5432,
            "postgresql_db": "reframe_benchmarks",
            "postgresql_conn_timeout": 60,
        }
    ],
    "modes": [
        {
            "name": "maintenance",
            "options": [
                "--exec-policy=async",
                "--reservation=maintenance",
                "--save-log-files",
                "--tag=acceptance",
                "--timestamp=%F_%H-%M-%S",
            ],
        },
    ],
    "general": [
        {
            "check_search_path": ["checks/"], 
            "check_search_recursive": True
        }
    ],
}


# ============================================================================
# Helper Functions
# ============================================================================

def get_partition_info(partition_name: str = "booster"):
    """
    Retrieve partition information from the configuration.
    
    Args:
        partition_name: Name of the partition (default: booster)
        
    Returns:
        Dictionary with partition information including:
        - num_nodes: Total available nodes
        - num_cpus_per_node: CPUs per node
        - num_gpus_per_node: GPUs per node
        - processor_arch: Processor architecture
        - gpu_model: GPU model (if applicable)
    """
    leonardo_system = site_configuration["systems"][0]
    
    for partition in leonardo_system["partitions"]:
        if partition["name"] == partition_name:
            processor = partition.get("processor", {})
            devices = partition.get("devices", [])
            
            info = {
                "partition_name": partition["name"],
                "description": partition["descr"],
                "scheduler": partition["scheduler"],
                "launcher": partition["launcher"],
                "num_cpus_per_node": processor.get("num_cpus", 0),
                "num_cpus_per_socket": processor.get("num_cpus_per_socket", 0),
                "num_sockets": processor.get("num_sockets", 0),
                "processor_arch": processor.get("arch", "unknown"),
                "num_gpus_per_node": devices[0].get("num_devices", 0) if devices else 0,
                "gpu_model": devices[0].get("model", "N/A") if devices else "N/A",
                "gpu_arch": devices[0].get("arch", "N/A") if devices else "N/A",
                "max_jobs": partition.get("max_jobs", 1),
                "available_environs": partition.get("environs", []),
            }
            
            return info
    
    return None


def validate_scaling_config(partition_name: str, num_nodes: int, procs_per_node: int):
    """
    Validate scaling configuration against partition limits.
    
    Args:
        partition_name: Name of the partition
        num_nodes: Number of nodes requested
        procs_per_node: Processes per node
        
    Returns:
        Tuple (is_valid, error_message)
    """
    partition_info = get_partition_info(partition_name)
    
    if not partition_info:
        return False, f"Partition '{partition_name}' not found in configuration"
    
    max_cpus = partition_info["num_cpus_per_node"]
    
    if procs_per_node > max_cpus:
        return False, (
            f"Requested {procs_per_node} processes per node, "
            f"but partition '{partition_name}' only has {max_cpus} CPUs per node"
        )
    
    return True, "Configuration valid"


# ============================================================================
# Example Usage
# ============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("LEONARDO System Configuration")
    print("=" * 80)
    
    # Display partition information
    for partition in ["login", "booster", "dcgp"]:
        print(f"\nPartition: {partition}")
        print("-" * 80)
        info = get_partition_info(partition)
        if info:
            for key, value in info.items():
                print(f"  {key:25s}: {value}")
    
    # Validate a scaling configuration
    print("\n" + "=" * 80)
    print("Validation Example")
    print("=" * 80)
    is_valid, message = validate_scaling_config("booster", num_nodes=4, procs_per_node=32)
    print(f"Valid: {is_valid}")
    print(f"Message: {message}")
    
    # Show available environments
    print("\n" + "=" * 80)
    print("Available Environments")
    print("=" * 80)
    for env in site_configuration["environments"]:
        print(f"\n  {env['name']}:")
        print(f"    Compilers: {env['cc']}, {env['cxx']}, {env['ftn']}")
        print(f"    Features: {', '.join(env['features'])}")
        if env['modules']:
            print(f"    Modules: {', '.join(env['modules'])}")
