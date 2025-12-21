"""
Slurm scheduler backend with full sbatch integration.

This module uses the CENTRALIZED topology detection to automatically infer
hardware configuration without hardcoded values.

GPU-aware configuration with automatic hardware detection:
- Automatic GPU/CPU detection via Slurm queries and environment
- Proper CPU allocation vs MPI task count
- Cores per GPU calculation
- GPU binding script generation

IMPORTANT: OpenMPI Syntax
=========================
The syntax `--map-by ppr:X:node:PE=Y` is NOT valid for OpenMPI.
The correct equivalent that achieves the same behavior is:

    mpirun -np N --map-by ppr:X:node --bind-to core --cpus-per-proc Y

This module generates ONLY correct, working syntax.

Author: HPC-ScaleTest Contributors
"""

import os
import subprocess
import logging
import re
import time
from pathlib import Path
from typing import List, Optional, Dict, Any

from core.abstracts import SchedulerInterface
from core.config import JobConfig, ResourceConfig
from core.types import JobStatus, LARGE_JOB_THRESHOLD

# Import the new job configuration generator
try:
    from core.job_config_generator import (
        PartitionDetector, PartitionTopology, JobConfiguration,
        generate_job_config, generate_gpu_binding_script as gen_bind_script,
        generate_slurm_job_script, get_partition_detector
    )
    HAS_JOB_CONFIG_GENERATOR = True
except ImportError:
    HAS_JOB_CONFIG_GENERATOR = False

# Import centralized topology detection - THIS IS THE SINGLE SOURCE OF TRUTH
try:
    from core.topology import (
        TopologyDetector, NodeTopology, MPIMapping,
        get_topology_detector, GPUVendor, DetectionContext
    )
    from core.mpi_command import (
        MPICommandGenerator, MPIDetector, MPIInfo,
        generate_gpu_binding_script
    )
    HAS_TOPOLOGY_MODULE = True
except ImportError:
    HAS_TOPOLOGY_MODULE = False

# Legacy imports for backward compatibility
try:
    from utils.slurm_detector import get_slurm_detector
    HAS_SLURM_DETECTOR = True
except ImportError:
    HAS_SLURM_DETECTOR = False

try:
    from utils.gpu_detection import GPUDetector
    from utils.partition_validator import PartitionValidator
    HAS_GPU_DETECTION = True
except ImportError:
    HAS_GPU_DETECTION = False

try:
    from utils.advanced_gpu_manager import AdvancedGPUManager
    HAS_ADVANCED_GPU = True
except ImportError:
    HAS_ADVANCED_GPU = False

logger = logging.getLogger(__name__)


class SlurmScheduler(SchedulerInterface):
    """
    Slurm workload manager backend with automatic topology detection.
    
    This scheduler uses the centralized topology detection from core/topology.py
    to automatically configure jobs without hardcoded assumptions.
    
    Detection Flow:
    1. If inside a Slurm job: Use SLURM_* environment variables
    2. If on login node: Query partition via sinfo/scontrol
    3. Fall back to system introspection if needed
    
    The topology module derives:
    - CPU cores per node
    - GPUs per node
    - MPI ranks per node (= GPUs for GPU jobs, = cores for CPU jobs)
    - CPU cores per rank
    """
    
    def __init__(self, options: Optional[Dict[str, Any]] = None):
        """
        Initialize the Slurm scheduler.
        
        Args:
            options: Optional configuration dictionary with keys:
                - partition: Default partition name
                - account: Default account name
                - qos: Default quality of service
                - time_limit: Default time limit
        """
        self.options = options or {}
        self._topology_detector = None
        self._mpi_generator = None
        self._partition_detector = None
    
    @property
    def topology_detector(self) -> 'TopologyDetector':
        """Get or create topology detector."""
        if self._topology_detector is None and HAS_TOPOLOGY_MODULE:
            self._topology_detector = get_topology_detector()
        return self._topology_detector
    
    def _detect_topology_and_mapping(
        self,
        job_config: JobConfig,
        resource_config: ResourceConfig,
    ) -> Dict[str, Any]:
        """
        Detect hardware topology and compute MPI mapping.
        
        This method uses partition-aware detection to automatically infer:
        - CPU cores per node
        - GPUs per node
        
        Then derives:
        - MPI ranks per node = GPUs per node (for GPU jobs)
        - CPU cores per rank = CPU cores / ranks
        
        Detection Methods (priority order):
        1. PartitionDetector (new, partition-specific)
        2. TopologyDetector (legacy, general)
        3. Legacy fallback methods
        
        Returns:
            Dictionary with keys:
                - topology: PartitionTopology or NodeTopology
                - job_config_obj: JobConfiguration (if available)
                - partition: str
                - gpus_per_node: int
                - tasks_per_node: int (MPI ranks)
                - cores_per_task: int
                - cpu_cores_per_node: int
                - is_gpu_job: bool
        """
        # Get partition (from config, environment, or use default)
        partition = resource_config.partition
        if not partition:
            partition = os.environ.get('HPC_SCALETEST_PARTITION')
        
        if not partition:
            logger.warning(
                "No partition specified. Set via resource_config.partition or "
                "HPC_SCALETEST_PARTITION environment variable."
            )
        
        # Check for user overrides
        user_ranks = getattr(resource_config, 'actual_mpi_tasks', None)
        user_cores = getattr(resource_config, 'cores_per_task', None)
        is_gpu_job = getattr(resource_config, 'hardware_type', '') == 'gpu'
        gpus_config = getattr(resource_config, 'gpus_per_node', 0)
        
        # =================================================================
        # Method 1: Use new PartitionDetector (partition-aware)
        # =================================================================
        if HAS_JOB_CONFIG_GENERATOR and partition:
            try:
                detector = get_partition_detector()
                part_topology = detector.detect(partition)
                
                # Override GPUs if user configured more than detected
                if gpus_config > part_topology.gpus_per_node:
                    part_topology = PartitionTopology(
                        partition=partition,
                        cpu_cores_per_node=part_topology.cpu_cores_per_node,
                        gpus_per_node=gpus_config,
                        detection_method=part_topology.detection_method + " + user_override"
                    )
                
                # Create JobConfiguration
                job_cfg = JobConfiguration(
                    num_nodes=job_config.num_nodes,
                    topology=part_topology
                )
                
                logger.info(f"Partition-aware topology detection:")
                logger.info(f"  Partition: {partition}")
                logger.info(f"  Detection method: {part_topology.detection_method}")
                logger.info(f"  CPU cores per node: {part_topology.cpu_cores_per_node}")
                logger.info(f"  GPUs per node: {part_topology.gpus_per_node}")
                logger.info(f"  Derived MPI mapping:")
                logger.info(f"    Ranks per node: {part_topology.ranks_per_node}")
                logger.info(f"    Cores per rank: {part_topology.cores_per_rank}")
                
                return {
                    'topology': part_topology,
                    'job_config_obj': job_cfg,
                    'partition': partition,
                    'gpus_per_node': part_topology.gpus_per_node,
                    'cpu_cores_per_node': part_topology.cpu_cores_per_node,
                    'tasks_per_node': part_topology.ranks_per_node,
                    'cores_per_task': part_topology.cores_per_rank,
                    'is_gpu_job': part_topology.gpus_per_node > 0,
                }
                
            except Exception as e:
                logger.warning(f"PartitionDetector failed: {e}")
                logger.info("Falling back to TopologyDetector")
        
        # =================================================================
        # Method 2: Use TopologyDetector (legacy)
        # =================================================================
        if HAS_TOPOLOGY_MODULE and self.topology_detector:
            try:
                # Detect topology for this partition
                topology = self.topology_detector.detect(partition)
                
                # Override GPUs if configured but not detected
                if gpus_config > 0 and topology.gpus == 0:
                    # User configured GPUs but none detected - trust user
                    topology = NodeTopology(
                        cpu_cores=topology.cpu_cores,
                        gpus=gpus_config,
                        gpu_vendor=GPUVendor.NVIDIA,  # Assume NVIDIA if not detected
                        partition=partition,
                        detection_method=topology.detection_method + " + user_override"
                    )
                
                # Compute MPI mapping
                mapping = self.topology_detector.compute_mpi_mapping(
                    topology=topology,
                    num_nodes=job_config.num_nodes,
                    user_ranks_per_node=user_ranks,
                    user_cores_per_rank=user_cores,
                )
                
                logger.info(f"Topology detection complete:")
                logger.info(f"  Method: {topology.detection_method}")
                logger.info(f"  CPU cores: {topology.cpu_cores}")
                logger.info(f"  GPUs: {topology.gpus}")
                logger.info(f"  MPI ranks/node: {mapping.ranks_per_node}")
                logger.info(f"  Cores/rank: {mapping.cores_per_rank}")
                
                return {
                    'topology': topology,
                    'mapping': mapping,
                    'partition': partition,
                    'gpus_per_node': topology.gpus,
                    'cpu_cores_per_node': topology.cpu_cores,  # ADDED: For full node allocation
                    'tasks_per_node': mapping.ranks_per_node,
                    'cores_per_task': mapping.cores_per_rank,
                    'is_gpu_job': topology.gpus > 0,
                }
                
            except Exception as e:
                logger.warning(f"Centralized topology detection failed: {e}")
                logger.warning("Falling back to legacy detection")
        
        # ===== LEGACY FALLBACK =====
        # This code path is only used if the topology module is not available
        # or if detection fails. It preserves backward compatibility.
        
        gpus_per_node = resource_config.gpus_per_node
        procs_per_node = resource_config.procs_per_node
        cores_per_task = 1
        tasks_per_node = procs_per_node
        
        # Use legacy GPU manager if available
        if HAS_ADVANCED_GPU and is_gpu_job:
            try:
                gpu_manager = AdvancedGPUManager()
                gpu_node_config = gpu_manager.detect_gpu_node_config(partition)
                
                if gpu_node_config:
                    gpus_per_node = gpu_node_config.gpus_per_node
                    cores_per_task = gpu_node_config.cores_per_gpu
                    procs_per_node = gpu_node_config.cpus_per_node
                    tasks_per_node = gpus_per_node  # 1 task per GPU
                    
                    logger.info(f"Legacy GPU detection: {gpus_per_node} GPUs, "
                              f"{cores_per_task} cores/GPU")
            except Exception as e:
                logger.warning(f"Legacy GPU detection failed: {e}")
        
        # Basic GPU detection fallback
        elif HAS_GPU_DETECTION and is_gpu_job:
            try:
                gpu_detector = GPUDetector()
                gpu_info = gpu_detector.detect(partition_name=partition)
                
                if gpu_info:
                    gpus_per_node = gpu_info.count_per_node
                    tasks_per_node = gpus_per_node
                    cores_per_task = procs_per_node // gpus_per_node if procs_per_node > 0 else 1
            except Exception as e:
                logger.warning(f"Basic GPU detection failed: {e}")
        
        # GPU job: set reasonable defaults if no detection worked
        if is_gpu_job and gpus_per_node > 0:
            if tasks_per_node != gpus_per_node:
                tasks_per_node = gpus_per_node
                cores_per_task = procs_per_node // gpus_per_node if procs_per_node > 0 else 1
        
        return {
            'topology': None,
            'mapping': None,
            'partition': partition,
            'gpus_per_node': gpus_per_node,
            'cpu_cores_per_node': procs_per_node,  # ADDED: For full node allocation
            'tasks_per_node': tasks_per_node,
            'cores_per_task': cores_per_task,
            'is_gpu_job': gpus_per_node > 0,
        }
    
    def generate_job_script(
        self,
        job_config: JobConfig,
        resource_config: ResourceConfig,
        command: List[str],
        env_setup: List[str]
    ) -> str:
        """
        Generate Slurm batch script with automatic topology detection.
        
        CRITICAL: This method correctly separates:
        - SLURM resource allocation (full node: --ntasks-per-node=CPUs)
        - MPI execution (actual tasks: -np with ppr mapping)
        
        The generated script includes:
        - Correct SLURM directives requesting FULL node resources
        - Proper MPI rank and CPU allocation via mpirun
        - GPU binding script for GPU jobs
        - Environment setup and module loading
        """
        
        # Detect topology and compute mapping
        config = self._detect_topology_and_mapping(job_config, resource_config)
        
        partition = config['partition']
        gpus_per_node = config['gpus_per_node']
        
        # Get CPU cores per node (for full allocation)
        cpu_cores_per_node = config.get('cpu_cores_per_node', resource_config.procs_per_node)
        
        # Actual MPI tasks (derived from topology)
        mpi_ranks_per_node = config['tasks_per_node']
        cores_per_rank = config['cores_per_task']
        is_gpu_job = config['is_gpu_job']
        
        # Store computed values back in resource_config for the launcher
        resource_config.actual_mpi_tasks = mpi_ranks_per_node
        resource_config.cores_per_task = cores_per_rank
        
        # Begin script generation
        script = """#!/bin/bash
# =============================================================================
# SLURM Job Script - Generated by HPC-ScaleTest
# =============================================================================
# Topology Detection Method: {}
#
# IMPORTANT: Resource Allocation vs MPI Execution
# ------------------------------------------------
# SLURM allocates FULL node resources for accounting:
#   --ntasks-per-node={} (all CPUs)
#   --gres=gpu:{} (all GPUs)
#
# MPI executes with derived task configuration:
#   {} MPI ranks per node (1 per GPU for GPU jobs)
#   {} CPU cores per rank
# =============================================================================

""".format(
            config.get('topology').detection_method if config.get('topology') else 'legacy',
            cpu_cores_per_node,
            gpus_per_node,
            mpi_ranks_per_node,
            cores_per_rank
        )
        
        # Slurm directives
        script += f"#SBATCH --nodes={job_config.num_nodes}\n"
        if partition:
            script += f"#SBATCH --partition={partition}\n"
        if resource_config.qos:
            script += f"#SBATCH --qos={resource_config.qos}\n"
        
        # CRITICAL FIX: Use FULL CPU count for --ntasks-per-node (resource allocation)
        # NOT the actual MPI ranks (which is what we run)
        script += f"#SBATCH --ntasks-per-node={cpu_cores_per_node}\n"
        
        # GPU allocation (full node)
        if is_gpu_job and gpus_per_node > 0:
            script += f"#SBATCH --gres=gpu:{gpus_per_node}\n"
        
        # Log the configuration
        logger.info(f"SLURM resource allocation:")
        logger.info(f"  --ntasks-per-node={cpu_cores_per_node} (FULL node allocation)")
        logger.info(f"  --gres=gpu:{gpus_per_node}")
        logger.info(f"MPI execution configuration:")
        logger.info(f"  MPI ranks per node: {mpi_ranks_per_node}")
        logger.info(f"  Cores per rank: {cores_per_rank}")
        logger.info(f"  Total MPI ranks: {job_config.num_nodes * mpi_ranks_per_node}")
            
        if resource_config.exclusive:
            script += "#SBATCH --exclusive\n"
        
        # Account - use environment variable if not set
        account = resource_config.account or os.environ.get('HPC_SCALETEST_ACCOUNT')
        if account:
            script += f"#SBATCH -A {account}\n"
        
        script += f"#SBATCH --time={resource_config.time_limit}\n"
        script += f"#SBATCH --job-name={job_config.job_id}\n"
        script += f"#SBATCH -o {job_config.working_dir}/job.out\n"
        script += f"#SBATCH -e {job_config.working_dir}/job.err\n"
        
        if resource_config.mail_user:
            script += f"#SBATCH --mail-type=ALL\n"
            script += f"#SBATCH --mail-user={resource_config.mail_user}\n"
        
        if resource_config.constraint:
            script += f"#SBATCH --constraint={resource_config.constraint}\n"
        
        if resource_config.reservation:
            script += f"#SBATCH --reservation={resource_config.reservation}\n"
        
        if resource_config.memory_per_node:
            script += f"#SBATCH --mem={resource_config.memory_per_node}\n"
        
        # Large job handling
        if job_config.num_nodes > LARGE_JOB_THRESHOLD:
            script += "# Note: Large job (>{} nodes)\n".format(LARGE_JOB_THRESHOLD)
        
        script += "\n# Environment setup\n"
        for cmd in env_setup:
            script += f"{cmd}\n"
        
        script += '\necho "Loading modules"\n'
        
        script += "\n# Change to working directory\n"
        if job_config.working_dir:
            script += f"cd {job_config.working_dir}\n"
        
        script += "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}\n"
        
        # =================================================================
        # BINARY PATH AND GPU BINDING - CRITICAL FIX
        # =================================================================
        # BINARY must point to the directory containing the compiled binary,
        # NOT ${PWD}. The executable must be the actual binary, not bind.sh.
        # =================================================================
        
        command_copy = command.copy() if command else []
        binary_set = False
        binary_name = None
        binary_dir = None
        
        # Find the executable in the command
        binary_idx = None
        for idx, arg in enumerate(command_copy):
            # Skip launcher and flags
            if arg.startswith('-') or arg in ['srun', 'mpirun', 'mpiexec', 'bind.sh', './bind.sh', 'bind_gpu', './bind_gpu']:
                continue
            if arg.startswith('$') and '/' not in arg:
                continue
            if arg.isdigit():
                continue
            if '=' in arg or ':' in arg:
                continue
            
            # This looks like an executable path
            if '/' in arg or (arg and arg[0].isalpha()):
                binary_idx = idx
                break
        
        # Extract and set binary path
        if binary_idx is not None:
            from pathlib import Path
            full_binary_path = command_copy[binary_idx]
            
            # Handle $BINARY/ paths
            if full_binary_path.startswith('$BINARY/'):
                # Already properly formatted
                binary_name = full_binary_path.replace('$BINARY/', '')
                script += f"\n# Binary location (from command)\n"
                script += f"# Executable: {binary_name}\n"
                binary_set = True
            elif '/' in full_binary_path:
                # Full or relative path - extract directory and name
                binary_path_obj = Path(full_binary_path)
                
                # CRITICAL: Use absolute path, not ${PWD}
                if binary_path_obj.is_absolute():
                    binary_dir = str(binary_path_obj.parent)
                else:
                    # Convert relative to absolute using resolve()
                    try:
                        binary_dir = str(binary_path_obj.resolve().parent)
                    except:
                        # If resolve fails, use the path as-is with working dir
                        binary_dir = str((Path(job_config.working_dir) / binary_path_obj).parent)
                
                binary_name = binary_path_obj.name
                
                script += f"\n# Binary location\n"
                script += f"export BINARY={binary_dir}\n"
                script += f"# Executable: {binary_name}\n"
                command_copy[binary_idx] = f"$BINARY/{binary_name}"
                binary_set = True
                
                logger.info(f"Job script binary configuration:")
                logger.info(f"  BINARY={binary_dir}")
                logger.info(f"  Executable={binary_name}")
            else:
                # Just executable name - try to find it in install_dir or working_dir
                binary_name = full_binary_path
                
                # Check if we have an install_dir from the test
                install_dir = getattr(job_config, 'install_dir', None)
                if install_dir:
                    binary_dir = str(install_dir)
                else:
                    binary_dir = str(job_config.working_dir)
                
                script += f"\n# Binary location\n"
                script += f"export BINARY={binary_dir}\n"
                script += f"# Executable: {binary_name}\n"
                command_copy[binary_idx] = f"$BINARY/{binary_name}"
                binary_set = True
                
                logger.info(f"Job script binary configuration:")
                logger.info(f"  BINARY={binary_dir}")
                logger.info(f"  Executable={binary_name}")
        
        # If no binary was set but we're in a working directory, use it
        if not binary_set and job_config.working_dir:
            script += f"\n# Binary location (default)\n"
            script += f"export BINARY={job_config.working_dir}\n"
        
        # =================================================================
        # GPU BINDING SCRIPT
        # =================================================================
        # For GPU jobs, generate bind.sh that uses OMPI_COMM_WORLD_LOCAL_RANK
        # to set CUDA_VISIBLE_DEVICES
        # =================================================================
        
        needs_gpu_binding = is_gpu_job and gpus_per_node > 0
        has_bind_script = any(arg in ['./bind.sh', 'bind.sh'] for arg in command_copy)
        
        if needs_gpu_binding:
            script += '''
# =============================================================================
# GPU Binding Script
# =============================================================================
# Uses OMPI_COMM_WORLD_LOCAL_RANK to bind each MPI rank to a unique GPU
# =============================================================================
cat > bind.sh << 'BIND_EOF'
#!/bin/bash
# Determine local rank from available environment variables
if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
    LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ -n "$MPI_LOCALRANKID" ]; then
    LOCAL_RANK=$MPI_LOCALRANKID
elif [ -n "$PMI_LOCAL_RANK" ]; then
    LOCAL_RANK=$PMI_LOCAL_RANK
elif [ -n "$SLURM_LOCALID" ]; then
    LOCAL_RANK=$SLURM_LOCALID
else
    LOCAL_RANK=0
fi
export CUDA_VISIBLE_DEVICES=$LOCAL_RANK
exec "$@"
BIND_EOF
chmod +x bind.sh
'''
            
            # Insert bind.sh wrapper if not already present
            if not has_bind_script and binary_idx is not None:
                command_copy.insert(binary_idx, './bind.sh')
            # Try to get binary name from resource config or use default
            if binary_name:
                script += f"# Executable: {binary_name}\n"
        
        script += "\ndate\n"
        
        # Main command with timing (full launch_cmd)
        script += 'start_time="$(date -u +%s.%N)"\n'
        script += "# Execute command\n"
        script += " ".join(command_copy) + "\n\n"
        script += 'end_time="$(date -u +%s.%N)"\n'
        script += 'elapsed="$(bc <<<"$end_time-$start_time")"\n'
        script += 'echo "Total of $elapsed seconds elapsed for process"\n'
        
        script += "date\n"
        
        return script
    
    def submit_job(self, script_path: Path, max_retries: int = 3) -> str:
        """Submit job using sbatch with verification and retry logic."""
        logger.info(f"Attempting to submit job script: {script_path}")
        
        # Verify script exists and is readable
        if not script_path.exists():
            raise FileNotFoundError(f"Job script not found: {script_path}")
        
        if not script_path.is_file():
            raise ValueError(f"Path is not a file: {script_path}")
        
        # Check if sbatch is available
        try:
            which_result = subprocess.run(["which", "sbatch"], capture_output=True, text=True)
            if which_result.returncode == 0:
                sbatch_path = which_result.stdout.strip()
                logger.debug(f"sbatch found at: {sbatch_path}")
            else:
                logger.error("❌ sbatch not found in PATH")
                logger.error("Are you on a login node with Slurm access?")
                raise RuntimeError("sbatch command not available")
        except Exception as e:
            logger.error(f"Could not verify sbatch availability: {e}")
            raise
        
        # Try submission with retries
        last_error = None
        for attempt in range(1, max_retries + 1):
            try:
                if attempt > 1:
                    wait_time = 2 ** (attempt - 1)  # Exponential backoff
                    logger.info(f"Retry attempt {attempt}/{max_retries} (waiting {wait_time}s)")
                    time.sleep(wait_time)
                
                logger.info(f"Executing: sbatch {script_path}")
                
                # Method 1: Direct sbatch call
                result = subprocess.run(
                    ["sbatch", str(script_path)],
                    capture_output=True,
                    text=True,
                    check=False  # Don't raise immediately, we'll check manually
                )
                
                # If direct call fails, try with shell wrapper (for environment loading)
                if result.returncode != 0:
                    logger.warning(f"Direct sbatch failed (exit {result.returncode}), trying with shell wrapper")
                    result = subprocess.run(
                        ["bash", "-lc", f"sbatch {script_path}"],
                        capture_output=True,
                        text=True,
                        check=False
                    )
                
                # Log all output for debugging
                logger.debug(f"sbatch return code: {result.returncode}")
                logger.debug(f"sbatch stdout: {result.stdout}")
                if result.stderr:
                    logger.debug(f"sbatch stderr: {result.stderr}")
                
                # Check return code
                if result.returncode != 0:
                    error_msg = f"sbatch failed with exit code {result.returncode}"
                    if result.stderr:
                        error_msg += f"\nstderr: {result.stderr}"
                    if result.stdout:
                        error_msg += f"\nstdout: {result.stdout}"
                    logger.error(error_msg)
                    last_error = RuntimeError(error_msg)
                    continue  # Try again
                
                # Parse job ID from stdout
                match = re.search(r'Submitted batch job (\d+)', result.stdout)
                if not match:
                    error_msg = f"Could not parse job ID from sbatch output: {result.stdout}"
                    logger.error(error_msg)
                    last_error = ValueError(error_msg)
                    continue  # Try again
                
                job_id = match.group(1)
                logger.info(f"✓ sbatch returned job ID: {job_id}")
                
                # Verify job was actually submitted by querying scheduler
                logger.debug(f"Verifying job {job_id} in queue...")
                verify_result = self._verify_job_submission(job_id)
                
                if verify_result:
                    logger.info(f"✓ Job {job_id} verified in Slurm queue")
                    return job_id
                else:
                    logger.warning(f"⚠ Job {job_id} not found in queue immediately (may be normal)")
                    # Still return the job ID - it might just be processing
                    return job_id
                    
            except subprocess.CalledProcessError as e:
                last_error = e
                logger.error(f"Attempt {attempt}/{max_retries} failed: {e}")
                logger.error(f"  Return code: {e.returncode}")
                logger.error(f"  stdout: {e.stdout}")
                logger.error(f"  stderr: {e.stderr}")
                
            except Exception as e:
                last_error = e
                logger.error(f"Attempt {attempt}/{max_retries} failed: {e}")
        
        # All retries exhausted
        logger.error(f"❌ Job submission failed after {max_retries} attempts")
        if last_error:
            raise last_error
        else:
            raise RuntimeError("Job submission failed for unknown reason")
    
    def _verify_job_submission(self, job_id: str) -> bool:
        """Verify job exists in Slurm queue."""
        try:
            # Check with squeue first
            result = subprocess.run(
                ["squeue", "-j", job_id, "-h"],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0 and result.stdout.strip():
                return True
            
            # Also check sacct (for very fast jobs)
            result = subprocess.run(
                ["sacct", "-j", job_id, "--noheader", "--format=JobID"],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            return result.returncode == 0 and job_id in result.stdout
            
        except Exception as e:
            logger.debug(f"Could not verify job {job_id}: {e}")
            return False

    def get_job_status(self, job_id: str) -> JobStatus:
        """Query job status using squeue/sacct."""
        try:
            # squeue first
            result = subprocess.run(
                ["squeue", "-j", job_id, "-h", "-o", "%T"],
                capture_output=True,
                text=True
            )
            
            if result.returncode == 0 and result.stdout.strip():
                state = result.stdout.strip().upper()
                return self._parse_slurm_state(state)
            
            # sacct for completed
            result = subprocess.run(
                ["sacct", "-j", job_id, "-n", "-o", "State"],
                capture_output=True,
                text=True
            )
            
            if result.returncode == 0 and result.stdout.strip():
                state = result.stdout.strip().split()[0].upper()
                return self._parse_slurm_state(state)
            
            return JobStatus.UNKNOWN
            
        except Exception as e:
            logger.error(f"Failed to query job status: {e}")
            return JobStatus.UNKNOWN
    
    def _parse_slurm_state(self, state: str) -> JobStatus:
        """Convert Slurm state to JobStatus."""
        state_map = {
            'PENDING': JobStatus.PENDING,
            'RUNNING': JobStatus.RUNNING,
            'COMPLETED': JobStatus.COMPLETED,
            'FAILED': JobStatus.FAILED,
            'CANCELLED': JobStatus.CANCELLED,
            'TIMEOUT': JobStatus.TIMEOUT,
            'NODE_FAIL': JobStatus.FAILED,
            'PREEMPTED': JobStatus.CANCELLED,
        }
        return state_map.get(state, JobStatus.UNKNOWN)
    
    def cancel_job(self, job_id: str) -> bool:
        """Cancel job using scancel."""
        try:
            subprocess.run(
                ["scancel", job_id],
                check=True,
                capture_output=True
            )
            logger.info(f"Cancelled job {job_id}")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"scancel failed: {e.stderr}")
            return False
    
    def wait_for_completion(
        self,
        job_id: str,
        timeout: Optional[int] = None
    ) -> JobStatus:
        """Wait for job to complete."""
        start_time = time.time()
        while True:
            status = self.get_job_status(job_id)
            if status not in (JobStatus.PENDING, JobStatus.RUNNING):
                return status
            if timeout and (time.time() - start_time) > timeout:
                return JobStatus.TIMEOUT
            time.sleep(10)
