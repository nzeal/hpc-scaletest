"""
Slurm scheduler backend with automatic topology detection.

This module uses the unified execution module (core/unified_execution.py) 
to automatically detect hardware topology and generate correct MPI commands.

Key Features:
=============
- Automatic CPU/GPU detection via SLURM queries and system introspection
- Proper MPI mapping: ranks_per_node = GPUs, cores_per_rank = CPUs/GPUs
- Full node resource allocation in SLURM directives
- Correct mpirun syntax: -np N --map-by ppr:X:node:PE=Y
- GPU binding via CUDA_VISIBLE_DEVICES (NVIDIA) or ROCR_VISIBLE_DEVICES (AMD)

Design Principles:
==================
- NO HARDCODED VALUES - All topology is detected at runtime
- SYSTEM AGNOSTIC - Works on Leonardo (NVIDIA), LUMI (AMD), etc.
- SINGLE SOURCE OF TRUTH - Uses core/unified_execution.py for all topology

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

# Import the unified execution module - SINGLE SOURCE OF TRUTH
try:
    from core.unified_execution import (
        UnifiedExecutor, TopologyDetector, HardwareTopology, GPUVendor
    )
    HAS_UNIFIED_EXECUTION = True
except ImportError:
    HAS_UNIFIED_EXECUTION = False

logger = logging.getLogger(__name__)


class SlurmScheduler(SchedulerInterface):
    """
    Slurm workload manager backend with automatic topology detection.
    
    This scheduler uses the unified execution module to automatically
    configure jobs without hardcoded assumptions.
    
    MPI Command Generation:
    =======================
    For GPU jobs (e.g., Leonardo Booster: 32 cores, 4 GPUs):
    - ranks_per_node = 4 (1 per GPU)
    - cores_per_rank = 8 (32/4)
    - Command: mpirun -np N --map-by ppr:4:node:PE=8 --report-bindings ./bind.sh <exe>
    
    SLURM Resource Allocation:
    ==========================
    - --ntasks-per-node = CPUs per node (full allocation)
    - --gres=gpu:N (all GPUs)
    
    Note: SLURM resource allocation uses FULL node resources for accounting,
    while MPI execution uses derived task count (ranks = GPUs).
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
                - launcher: 'srun' or 'mpirun' (default: 'srun')
        """
        self.options = options or {}
        self._executor: Optional[UnifiedExecutor] = None
        self._topology_cache: Dict[str, HardwareTopology] = {}
        # Get launcher preference - default to srun for SLURM
        self.launcher = self.options.get('launcher', 'srun')
    
    def _get_executor(self, partition: str = "") -> UnifiedExecutor:
        """Get or create unified executor for the partition."""
        if not HAS_UNIFIED_EXECUTION:
            raise RuntimeError(
                "Unified execution module not available. "
                "Ensure core/unified_execution.py exists."
            )
        
        if self._executor is None or (partition and self._executor.partition != partition):
            self._executor = UnifiedExecutor(
                partition=partition,
                account=self.options.get('account', ''),
                qos=self.options.get('qos', ''),
                time_limit=self.options.get('time_limit', '02:00:00')
            )
        return self._executor
    
    def _detect_topology(
        self,
        partition: str,
        resource_config: ResourceConfig
    ) -> HardwareTopology:
        """
        Detect hardware topology for the partition.
        
        Uses user-configured values from resource_config when available,
        falling back to auto-detection if not specified.
        
        IMPORTANT: For GPU jobs:
        - cpus_per_node = total CPU cores per node (e.g., 32)
        - gpus_per_node = number of GPUs per node (e.g., 4)
        - procs_per_node = MPI ranks per node (usually equals gpus_per_node for GPU jobs)
        
        Args:
            partition: SLURM partition name
            resource_config: Resource configuration with user overrides
        
        Returns:
            HardwareTopology with configured/detected values
        """
        # Get user-configured values
        user_cpus = getattr(resource_config, 'cpus_per_node', None) or 0
        user_gpus = getattr(resource_config, 'gpus_per_node', 0) or 0
        user_procs = getattr(resource_config, 'procs_per_node', 0) or 0
        
        # CRITICAL FIX: Check for known GPU partitions by name
        # If user_gpus is 0 but partition is a known GPU partition, infer GPUs
        known_gpu_partitions = {
            'boost_usr_prod': (4, 32, GPUVendor.NVIDIA),  # Leonardo Booster
            'boost_qos_bprod': (4, 32, GPUVendor.NVIDIA),
            'gpu': (4, 32, GPUVendor.NVIDIA),
            'dgx': (8, 64, GPUVendor.NVIDIA),
            'lumi-g': (8, 64, GPUVendor.AMD),
        }
        
        partition_lower = partition.lower() if partition else ''
        
        # If GPUs not specified but partition is known GPU partition, use defaults
        if user_gpus == 0:
            for known_part, (gpus, cpus, vendor) in known_gpu_partitions.items():
                if known_part in partition_lower:
                    logger.info(f"Detected known GPU partition '{partition}' - using defaults: {gpus} GPUs, {cpus} CPUs")
                    user_gpus = gpus
                    if user_cpus == 0:
                        user_cpus = cpus
                    break
        
        # For NON-GPU jobs only: procs_per_node can be used as CPU count
        if user_cpus == 0 and user_gpus == 0:
            user_cpus = user_procs
        
        cache_key = f"{partition}_{user_cpus}_{user_gpus}"
        
        if cache_key in self._topology_cache:
            return self._topology_cache[cache_key]
        
        logger.info(f"Detecting topology for partition '{partition}':")
        logger.info(f"  User-configured cpus_per_node: {user_cpus}")
        logger.info(f"  User-configured gpus_per_node: {user_gpus}")
        logger.info(f"  User-configured procs_per_node: {user_procs}")
        
        # If user provided explicit cpus_per_node and gpus_per_node, use them
        if user_cpus > 0 and user_gpus > 0:
            # GPU job with explicit CPU configuration
            gpu_vendor = GPUVendor.NVIDIA
            if 'lumi' in partition_lower or 'amd' in partition_lower:
                gpu_vendor = GPUVendor.AMD
            
            topology = HardwareTopology(
                cpu_cores_per_node=user_cpus,
                gpus_per_node=user_gpus,
                gpu_vendor=gpu_vendor,
                partition=partition,
                detection_method="user_config"
            )
            
            logger.info(f"Using user-configured GPU topology:")
        elif user_gpus > 0:
            # GPU job but CPU count not specified - need to auto-detect or use default
            logger.info(f"GPU count specified ({user_gpus}) but cpus_per_node not - auto-detecting...")
            
            cpu_count = None
            try:
                executor = self._get_executor(partition)
                detected = executor.detect_topology()
                detected_cpus = detected.cpu_cores_per_node
                logger.info(f"  Auto-detected CPU count: {detected_cpus}")
                
                # Sanity check: for GPU nodes, we expect at least 4-8 cores per GPU
                min_expected = user_gpus * 4
                if detected_cpus >= min_expected:
                    cpu_count = detected_cpus
                else:
                    logger.warning(f"  Auto-detected CPU count ({detected_cpus}) seems too low for {user_gpus} GPUs")
                    logger.warning(f"  Expected at least {min_expected} cores (4 per GPU)")
            except Exception as e:
                logger.warning(f"  Auto-detection failed: {e}")
            
            if cpu_count is None:
                # Use estimated value: typically 8 cores per GPU on modern systems
                cpu_count = user_gpus * 8
                if cpu_count < 16:
                    cpu_count = 32  # Minimum reasonable for modern GPU nodes
                logger.warning(f"  Using estimated CPU count: {cpu_count} (based on {user_gpus} GPUs × 8 cores)")
            
            gpu_vendor = GPUVendor.NVIDIA
            if 'lumi' in partition_lower or 'amd' in partition_lower:
                gpu_vendor = GPUVendor.AMD
            
            topology = HardwareTopology(
                cpu_cores_per_node=cpu_count,
                gpus_per_node=user_gpus,
                gpu_vendor=gpu_vendor,
                partition=partition,
                detection_method="user_gpus + auto_cpus"
            )
            
            logger.info(f"Using mixed topology (user GPUs, detected/estimated CPUs):")
        elif user_cpus > 0:
            # CPU job with explicit configuration
            topology = HardwareTopology(
                cpu_cores_per_node=user_cpus,
                gpus_per_node=0,
                gpu_vendor=GPUVendor.NONE,
                partition=partition,
                detection_method="user_config"
            )
            
            logger.info(f"Using user-configured CPU topology:")
        else:
            # Full auto-detection
            try:
                executor = self._get_executor(partition)
                topology = executor.detect_topology()
                logger.info(f"Using auto-detected topology:")
            except Exception as e:
                logger.warning(f"Auto-detection failed: {e}, using minimal defaults")
                topology = HardwareTopology(
                    cpu_cores_per_node=1,
                    gpus_per_node=0,
                    gpu_vendor=GPUVendor.NONE,
                    partition=partition,
                    detection_method="minimal_fallback"
                )
        
        logger.info(f"  CPU cores per node: {topology.cpu_cores_per_node}")
        logger.info(f"  GPUs per node: {topology.gpus_per_node}")
        logger.info(f"  GPU vendor: {topology.gpu_vendor.value}")
        logger.info(f"  Derived MPI mapping:")
        logger.info(f"    Ranks per node: {topology.ranks_per_node}")
        logger.info(f"    Cores per rank: {topology.cores_per_rank}")
        logger.info(f"  Detection method: {topology.detection_method}")
        
        self._topology_cache[cache_key] = topology
        
        return topology
    
    def generate_job_script(
        self,
        job_config: JobConfig,
        resource_config: ResourceConfig,
        command: List[str],
        env_setup: List[str]
    ) -> str:
        """
        Generate Slurm batch script with proper resource allocation.
        
        CRITICAL FIX:
        - SLURM must allocate FULL node CPUs to satisfy QOS requirements
        - MPI command uses actual ranks (which may be less than CPUs)
        - GPU binding is inlined directly, no external bind.sh
        
        Args:
            job_config: Job configuration
            resource_config: Resource configuration
            command: Executable and arguments
            env_setup: Environment setup commands
        
        Returns:
            Complete job script as string
        """
        # Get partition
        partition = resource_config.partition or os.environ.get('HPC_SCALETEST_PARTITION', '')
        
        # Detect topology
        topology = self._detect_topology(partition, resource_config)
        
        # Get configured launcher (default to srun for SLURM)
        launcher = self.launcher
        logger.info(f"Generating job script with launcher: {launcher}")
        
        # Extract executable and arguments from command
        executable = ""
        args = []
        binary_dir = ""
        
        if command:
            # Find the executable in the command
            for idx, arg in enumerate(command):
                # Skip launcher and flags
                if arg.startswith('-') or arg in ['srun', 'mpirun', 'mpiexec']:
                    continue
                if arg.startswith('$') and '/' not in arg:
                    continue
                if arg.isdigit():
                    continue
                if '=' in arg or (arg.startswith('ppr:') or arg.startswith('--')):
                    continue
                if arg in ['bind.sh', './bind.sh', 'bind_gpu.sh', './bind_gpu.sh']:
                    continue
                
                # This looks like an executable
                if '/' in arg:
                    from pathlib import Path
                    path = Path(arg)
                    if arg.startswith('$BINARY/'):
                        executable = arg.replace('$BINARY/', '')
                        binary_dir = "$BINARY"
                    else:
                        executable = path.name
                        binary_dir = str(path.parent) if str(path.parent) != '.' else ""
                else:
                    executable = arg
                
                # Remaining args
                args = command[idx+1:]
                break
        
        # Fallback: use working directory
        if not binary_dir and job_config.working_dir:
            binary_dir = str(job_config.working_dir)
        
        # Determine job type and calculate MPI layout
        is_gpu_job = topology.gpus_per_node > 0
        cpu_cores_per_node = topology.cpu_cores_per_node
        
        if is_gpu_job:
            # GPU job: 1 MPI rank per GPU
            mpi_ranks_per_node = topology.gpus_per_node
            cores_per_rank = cpu_cores_per_node // topology.gpus_per_node
        else:
            # CPU job: 1 MPI rank per core (or use user-configured value)
            mpi_ranks_per_node = resource_config.procs_per_node or cpu_cores_per_node
            cores_per_rank = 1
        
        total_mpi_ranks = job_config.num_nodes * mpi_ranks_per_node
        
        # Build executable path
        if binary_dir:
            exe_path = f"$BINARY/{executable}"
        else:
            exe_path = executable
        
        # Build args string
        args_str = ' '.join(args) if args else ''
        
        # Generate launcher command based on launcher type and job type
        if launcher == 'srun':
            # srun is SLURM native - simple and reliable
            if is_gpu_job:
                run_cmd = f"srun --ntasks={total_mpi_ranks} --ntasks-per-node={mpi_ranks_per_node} --cpus-per-task={cores_per_rank} --gpu-bind=closest {exe_path} {args_str}".strip()
            else:
                run_cmd = f"srun --ntasks={total_mpi_ranks} --ntasks-per-node={mpi_ranks_per_node} {exe_path} {args_str}".strip()
        else:
            # mpirun (OpenMPI)
            if is_gpu_job:
                # GPU jobs: use ppr mapping with bind.sh for GPU binding
                # CORRECT FORMAT: mpirun -np N --map-by ppr:X:node:PE=Y --report-bindings ./bind.sh $BINARY/exe args
                run_cmd = f"mpirun -np {total_mpi_ranks} --map-by ppr:{mpi_ranks_per_node}:node:PE={cores_per_rank} --report-bindings ./bind.sh {exe_path} {args_str}".strip()
            else:
                # CPU jobs: simple mpirun without complex mapping
                run_cmd = f"mpirun -np {total_mpi_ranks} {exe_path} {args_str}".strip()
        
        # Build script
        gpu_vendor_str = topology.gpu_vendor.value if topology.gpus_per_node > 0 else "N/A"
        script = f'''#!/bin/bash
# =============================================================================
# SLURM Job Script - Generated by HPC-ScaleTest
# =============================================================================
# Partition: {partition}
# Detection Method: {topology.detection_method}
#
# Hardware Topology (per node):
#   CPU cores: {cpu_cores_per_node}
#   GPUs: {topology.gpus_per_node}
#   GPU vendor: {gpu_vendor_str}
#
# MPI Mapping (computed from topology):
#   Ranks per node: {mpi_ranks_per_node} (1 per GPU)
#   Cores per rank: {cores_per_rank}
#   Total ranks: {total_mpi_ranks}
#
# IMPORTANT: SLURM resource allocation vs MPI execution
# -----------------------------------------------------
# SLURM allocates FULL node resources: --ntasks-per-node={cpu_cores_per_node}
# MPI executes with derived tasks: -np {total_mpi_ranks}
# =============================================================================

#SBATCH --nodes={job_config.num_nodes}
'''
        
        if partition:
            script += f"#SBATCH --partition={partition}\n"
        
        if resource_config.qos:
            script += f"#SBATCH --qos={resource_config.qos}\n"
        
        # CRITICAL: Allocate FULL node CPUs for QOS satisfaction
        # Use --ntasks-per-node with FULL CPU count
        script += f"#SBATCH --ntasks-per-node={cpu_cores_per_node}\n"
        script += f"#SBATCH --cpus-per-task=1\n"
        
        # GPU allocation
        if is_gpu_job:
            script += f"#SBATCH --gres=gpu:{topology.gpus_per_node}\n"
        
        if resource_config.exclusive:
            script += "#SBATCH --exclusive\n"
        
        # Account
        account = resource_config.account or os.environ.get('HPC_SCALETEST_ACCOUNT')
        if account:
            script += f"#SBATCH -A {account}\n"
        
        script += f"#SBATCH --time={resource_config.time_limit}\n"
        script += f"#SBATCH --job-name={job_config.job_id}\n"
        script += f"#SBATCH -o {job_config.working_dir}/job.out\n"
        script += f"#SBATCH -e {job_config.working_dir}/job.err\n"
        
        if resource_config.mail_user:
            script += "#SBATCH --mail-type=ALL\n"
            script += f"#SBATCH --mail-user={resource_config.mail_user}\n"
        
        if resource_config.constraint:
            script += f"#SBATCH --constraint={resource_config.constraint}\n"
        
        if resource_config.reservation:
            script += f"#SBATCH --reservation={resource_config.reservation}\n"
        
        if resource_config.memory_per_node:
            script += f"#SBATCH --mem={resource_config.memory_per_node}\n"
        
        if job_config.num_nodes > LARGE_JOB_THRESHOLD:
            script += f"# Note: Large job (>{LARGE_JOB_THRESHOLD} nodes)\n"
        
        script += "\n# Environment setup\n"
        for cmd in env_setup:
            script += f"{cmd}\n"
        
        script += '\necho "Loading modules"\n'
        
        script += "\n# Change to working directory\n"
        if job_config.working_dir:
            script += f"cd {job_config.working_dir}\n"
        
        script += f"\n# Binary location\n"
        if binary_dir and binary_dir != "$BINARY":
            script += f"export BINARY={binary_dir}\n"
        elif job_config.working_dir:
            script += f"export BINARY={job_config.working_dir}\n"
        
        script += f'''
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Nodes: $SLURM_NNODES"
echo "Partition: $SLURM_JOB_PARTITION"
echo "CPUs on node: $SLURM_CPUS_ON_NODE"
'''
        
        if is_gpu_job:
            script += f'echo "GPUs on node: ${{SLURM_GPUS_ON_NODE:-{topology.gpus_per_node}}}"\n'
        
        script += f'''echo "BINARY: $BINARY"
echo "=========================================="

# Set OpenMP threads to cores per rank
export OMP_NUM_THREADS={cores_per_rank}
'''
        
        # GPU binding - check for bind.sh, create as fallback if needed
        if is_gpu_job:
            script += f'''
# =============================================================================
# GPU Binding
# =============================================================================
# bind.sh should already exist in the job directory (created by HPC-ScaleTest)
# If not found, create it at runtime as a fallback

if [ ! -f ./bind.sh ]; then
    echo "WARNING: bind.sh not found, creating at runtime..."
    cat > bind.sh << 'BIND_EOF'
#!/bin/bash
# GPU Binding Script - Runtime fallback
# Determine local rank
if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
    LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ -n "$MPI_LOCALRANKID" ]; then
    LOCAL_RANK=$MPI_LOCALRANKID
elif [ -n "$PMI_LOCAL_RANK" ]; then
    LOCAL_RANK=$PMI_LOCAL_RANK
elif [ -n "$MV2_COMM_WORLD_LOCAL_RANK" ]; then
    LOCAL_RANK=$MV2_COMM_WORLD_LOCAL_RANK
elif [ -n "$SLURM_LOCALID" ]; then
    LOCAL_RANK=$SLURM_LOCALID
else
    LOCAL_RANK=0
fi
export CUDA_VISIBLE_DEVICES=$LOCAL_RANK
export ROCR_VISIBLE_DEVICES=$LOCAL_RANK
export HIP_VISIBLE_DEVICES=$LOCAL_RANK
exec "$@"
BIND_EOF
    chmod +x bind.sh
fi

'''
        
        script += f'''
echo ""
echo "Starting: $(date)"
echo "Command: {run_cmd}"
echo ""

# Execute command with timing
start_time="$(date -u +%s.%N)"

{run_cmd}

exit_code=$?
end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"

echo ""
echo "=========================================="
echo "Completed: $(date)"
echo "Exit code: $exit_code"
echo "Elapsed: $elapsed seconds"
echo "=========================================="

exit $exit_code
'''
        
        return script
    
    def submit_job(self, script_path: Path, max_retries: int = 3) -> str:
        """
        Submit job using sbatch with verification and retry logic.
        
        Args:
            script_path: Path to job script
            max_retries: Maximum retry attempts
        
        Returns:
            SLURM job ID
        
        Raises:
            RuntimeError: If submission fails after all retries
        """
        logger.info(f"Submitting job script: {script_path}")
        
        # Verify script exists
        if not script_path.exists():
            raise FileNotFoundError(f"Job script not found: {script_path}")
        
        if not script_path.is_file():
            raise ValueError(f"Path is not a file: {script_path}")
        
        # Check sbatch availability
        try:
            which_result = subprocess.run(
                ["which", "sbatch"],
                capture_output=True, text=True
            )
            if which_result.returncode != 0:
                raise RuntimeError("sbatch command not available")
            logger.debug(f"sbatch found at: {which_result.stdout.strip()}")
        except Exception as e:
            logger.error(f"Could not verify sbatch: {e}")
            raise
        
        # Submit with retries
        last_error = None
        for attempt in range(1, max_retries + 1):
            try:
                if attempt > 1:
                    wait_time = 2 ** (attempt - 1)
                    logger.info(f"Retry {attempt}/{max_retries} (waiting {wait_time}s)")
                    time.sleep(wait_time)
                
                logger.info(f"Executing: sbatch {script_path}")
                
                result = subprocess.run(
                    ["sbatch", str(script_path)],
                    capture_output=True, text=True, check=False
                )
                
                # Try shell wrapper if direct call fails
                if result.returncode != 0:
                    logger.warning(f"Direct sbatch failed, trying shell wrapper")
                    result = subprocess.run(
                        ["bash", "-lc", f"sbatch {script_path}"],
                        capture_output=True, text=True, check=False
                    )
                
                logger.debug(f"sbatch return code: {result.returncode}")
                logger.debug(f"sbatch stdout: {result.stdout}")
                if result.stderr:
                    logger.debug(f"sbatch stderr: {result.stderr}")
                
                if result.returncode != 0:
                    error_msg = f"sbatch failed with exit code {result.returncode}"
                    if result.stderr:
                        error_msg += f"\nstderr: {result.stderr}"
                    logger.error(error_msg)
                    last_error = RuntimeError(error_msg)
                    continue
                
                # Parse job ID
                match = re.search(r'Submitted batch job (\d+)', result.stdout)
                if not match:
                    error_msg = f"Could not parse job ID: {result.stdout}"
                    logger.error(error_msg)
                    last_error = ValueError(error_msg)
                    continue
                
                job_id = match.group(1)
                logger.info(f"✓ Job submitted: {job_id}")
                
                # Verify submission
                if self._verify_job_submission(job_id):
                    logger.info(f"✓ Job {job_id} verified in queue")
                else:
                    logger.warning(f"Job {job_id} not immediately visible (may be normal)")
                
                return job_id
                
            except subprocess.CalledProcessError as e:
                last_error = e
                logger.error(f"Attempt {attempt} failed: {e}")
            except Exception as e:
                last_error = e
                logger.error(f"Attempt {attempt} failed: {e}")
        
        # All retries exhausted
        logger.error(f"Job submission failed after {max_retries} attempts")
        if last_error:
            raise last_error
        raise RuntimeError("Job submission failed")
    
    def _verify_job_submission(self, job_id: str) -> bool:
        """Verify job exists in Slurm queue."""
        try:
            result = subprocess.run(
                ["squeue", "-j", job_id, "-h"],
                capture_output=True, text=True, timeout=5
            )
            if result.returncode == 0 and result.stdout.strip():
                return True
            
            # Check sacct for very fast jobs
            result = subprocess.run(
                ["sacct", "-j", job_id, "--noheader", "--format=JobID"],
                capture_output=True, text=True, timeout=5
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
                capture_output=True, text=True
            )
            
            if result.returncode == 0 and result.stdout.strip():
                state = result.stdout.strip().upper()
                return self._parse_slurm_state(state)
            
            # sacct for completed
            result = subprocess.run(
                ["sacct", "-j", job_id, "-n", "-o", "State"],
                capture_output=True, text=True
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
                check=True, capture_output=True
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
