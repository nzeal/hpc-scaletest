"""
GPU Scaling Manager - HARDENED VERSION

Manages GPU-aware scaling studies, supporting multi-vendor GPUs (NVIDIA, AMD, Intel)
and various MPI implementations (OpenMPI, MPICH, Cray).

CRITICAL FIXES:
- Removed hardcoded assumptions about NVIDIA/CUDA.
- Implemented dynamic environment variable generation for GPU visibility.
- Added logic for different MPI/GPU binding mechanisms.
"""

import logging
from typing import Dict, List, Optional, Tuple, Any

logger = logging.getLogger(__name__)

# --- Configuration for GPU/MPI Mapping ---
# This configuration centralizes the logic for different vendor/MPI combinations.
# It would ideally be loaded from a system config file.

GPU_VENDOR_MAP = {
    "nvidia": {
        "env_var": "CUDA_VISIBLE_DEVICES",
        "mpi_bind_flags": {
            "openmpi": "--bind-to none",  # OpenMPI often requires explicit unbinding
            "mpich": "--bind-to none",    # MPICH often requires explicit unbinding
            "cray": "--gpu-bind=closest", # Cray MPI uses specific flags
            "intel": "--bind-to none",
        },
        "srun_bind_flag": "--gpu-bind=closest",
    },
    "amd": {
        "env_var": "ROCR_VISIBLE_DEVICES",
        "mpi_bind_flags": {
            "openmpi": "--bind-to none",
            "mpich": "--bind-to none",
            "cray": "--gpu-bind=closest",
            "intel": "--bind-to none",
        },
        "srun_bind_flag": "--gpu-bind=closest",
    },
    "intel": {
        "env_var": "ZE_AFFINITY_MASK", # Or similar for oneAPI/Level Zero
        "mpi_bind_flags": {
            "openmpi": "--bind-to none",
            "mpich": "--bind-to none",
            "cray": "--gpu-bind=closest",
            "intel": "--bind-to none",
        },
        "srun_bind_flag": "--gpu-bind=closest",
    },
}

class GPUScalingManager:
    """
    Manages the environment and command line arguments for GPU-aware scaling.
    """

    def __init__(self, system_config: Dict[str, Any]):
        """
        Initializes the manager with the system configuration.
        
        Args:
            system_config: A dictionary containing system details, including
                           GPU vendor and MPI implementation.
        """
        self.system_config = system_config
        self.gpu_vendor = self._get_gpu_vendor()
        self.mpi_implementation = self._get_mpi_implementation()
        self.vendor_config = GPU_VENDOR_MAP.get(self.gpu_vendor, {})
        
        if not self.vendor_config:
            logger.warning(f"Unsupported GPU vendor: {self.gpu_vendor}. GPU scaling may not be optimal.")

    def _get_gpu_vendor(self) -> str:
        """Determines the GPU vendor from the system configuration."""
        # Assuming system_config has a 'devices' key with vendor info
        devices = self.system_config.get('devices', [])
        if devices:
            # Simple heuristic: take the vendor of the first device
            vendor = devices[0].get('vendor', 'unknown').lower()
            if 'nvidia' in vendor:
                return 'nvidia'
            elif 'amd' in vendor or 'rocm' in vendor:
                return 'amd'
            elif 'intel' in vendor or 'oneapi' in vendor:
                return 'intel'
        return 'unknown'

    def _get_mpi_implementation(self) -> str:
        """Determines the MPI implementation from the system configuration."""
        # Assuming system_config has a 'mpi' key
        mpi_impl = self.system_config.get('mpi', {}).get('implementation', 'unknown').lower()
        if 'openmpi' in mpi_impl:
            return 'openmpi'
        elif 'mpich' in mpi_impl and 'cray' in mpi_impl:
            return 'cray' # Cray MPICH is often handled differently
        elif 'mpich' in mpi_impl:
            return 'mpich'
        elif 'intel' in mpi_impl:
            return 'intel'
        return 'unknown'

    def get_gpu_environment_vars(self, local_rank: int, gpus_per_node: int) -> Dict[str, str]:
        """
        Generates environment variables to set GPU visibility for a specific process.
        
        Args:
            local_rank: The rank of the process within the node (0 to gpus_per_node - 1).
            gpus_per_node: The total number of GPUs available on the node.
            
        Returns:
            A dictionary of environment variables (e.g., {'CUDA_VISIBLE_DEVICES': '0'}).
        """
        env_vars = {}
        
        if not self.vendor_config:
            return env_vars

        env_var_name = self.vendor_config.get("env_var")
        
        if env_var_name:
            # The local_rank directly maps to the device index for most systems
            # (assuming the system has already set up the correct device ordering)
            device_index = str(local_rank)
            
            # For simplicity and robustness, we set the visible device to the local rank
            # This relies on the launcher (srun/mpirun) to handle the correct binding.
            env_vars[env_var_name] = device_index
            logger.debug(f"Setting {env_var_name}={device_index} for local rank {local_rank}")
            
        return env_vars

    def get_mpi_binding_flags(self) -> List[str]:
        """
        Generates MPI-specific binding flags to ensure correct process-to-GPU affinity.
        
        Returns:
            A list of command-line arguments for the MPI launcher (e.g., mpirun).
        """
        if not self.vendor_config:
            return []

        mpi_flags = self.vendor_config.get("mpi_bind_flags", {}).get(self.mpi_implementation)
        
        if mpi_flags:
            logger.debug(f"Using MPI binding flags for {self.mpi_implementation}: {mpi_flags}")
            return mpi_flags.split()
        
        logger.warning(f"No specific MPI binding flags found for {self.mpi_implementation} with {self.gpu_vendor}. Using default.")
        return []

    def get_srun_binding_flags(self) -> List[str]:
        """
        Generates Slurm-specific binding flags for GPU affinity.
        
        Returns:
            A list of command-line arguments for the srun launcher.
        """
        if not self.vendor_config:
            return []

        srun_flag = self.vendor_config.get("srun_bind_flag")
        
        if srun_flag:
            logger.debug(f"Using srun binding flag: {srun_flag}")
            return [srun_flag]
        
        return []

    def get_launcher_prefix(self, launcher: str) -> List[str]:
        """
        Generates the full command prefix (MPI flags + srun flags) for the launcher.
        
        Args:
            launcher: The name of the launcher ('mpirun', 'srun', etc.).
            
        Returns:
            A list of command-line arguments to prefix the application command.
        """
        prefix = []
        
        if launcher == 'srun':
            # Slurm handles both MPI and GPU binding via srun flags
            prefix.extend(self.get_srun_binding_flags())
        elif launcher == 'mpirun':
            # mpirun requires MPI-specific binding flags
            prefix.extend(self.get_mpi_binding_flags())
        
        # Note: The environment variables (CUDA_VISIBLE_DEVICES) are typically
        # set by the job script or the launcher itself, not passed as a prefix.
        # This function only returns command-line flags.
        
        return prefix

# Re-export the class for easy import
GPUScalingManager.GPU_VENDOR_MAP = GPU_VENDOR_MAP
