"""
Slurm scheduler backend with full sbatch integration.
Customized for user template: gres, exclusive, qos, mpirun/srun, timing, modules, etc.
"""

import subprocess
import logging
import re
import time
from pathlib import Path
from typing import List, Optional

from core.abstracts import SchedulerInterface
from core.config import JobConfig, ResourceConfig
from core.types import JobStatus, LARGE_JOB_THRESHOLD


logger = logging.getLogger(__name__)


class SlurmScheduler(SchedulerInterface):
    """Slurm workload manager backend."""
    
    def generate_job_script(
        self,
        job_config: JobConfig,
        resource_config: ResourceConfig,
        command: List[str],
        env_setup: List[str]
    ) -> str:
        """Generate Slurm batch script matching user template."""
        script = """#!/bin/bash

echo "======================================"
echo "======================================"
echo "This is my script"
cat $0
echo "======================================"
echo "======================================"

"""
        
        # Slurm directives
        script += f"#SBATCH --nodes={job_config.num_nodes}\n"
        script += f"#SBATCH --partition={resource_config.partition or 'X_usr_prod'}\n"
        if resource_config.qos:
            script += f"###SBATCH --qos={resource_config.qos}\n"
        script += f"#SBATCH --ntasks-per-node={resource_config.procs_per_node}\n"
        script += "#SBATCH --cpus-per-task=1\n"
        if resource_config.gpus_per_node > 0:
            script += f"#SBATCH --gres=gpu:{resource_config.gpus_per_node}\n"
        if resource_config.exclusive:
            script += "###SBATCH --exclusive\n"
        script += f"#SBATCH -A {resource_config.account or 'cin_X'}\n"
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
            script += "# Large job: adjust qos if needed\n"
        
        script += "\n# Environment setup\n"
        for cmd in env_setup:
            script += f"{cmd}\n"
        
        script += '\necho "Loading modules"\n'
        
        script += "\n# Change to working directory\n"
        if job_config.working_dir:
            script += f"cd {job_config.working_dir}\n"
        
        script += "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}\n"
        
        # Extract binary path as variable for clarity
        # Make a copy to avoid modifying the original command list
        command_copy = command.copy() if command else []
        
        if command_copy and len(command_copy) > 0:
            # Find the executable in the command (must contain '/' and not be a flag)
            binary_idx = None
            for idx, arg in enumerate(command_copy):
                # Skip flags (anything starting with -)
                if arg.startswith('-'):
                    continue
                # Skip launcher commands
                if arg in ['srun', 'mpirun', 'mpiexec']:
                    continue
                # Must be a path (contains /) and not a number
                if '/' in arg and not arg.replace('.', '').replace('/', '').isdigit():
                    binary_idx = idx
                    break
            
            # If we found a binary path, extract directory and name as variables
            if binary_idx is not None:
                from pathlib import Path
                full_binary_path = command_copy[binary_idx]
                binary_path_obj = Path(full_binary_path)
                binary_dir = str(binary_path_obj.parent)
                binary_name = binary_path_obj.name
                
                script += f"\n# Binary location\n"
                script += f"BINARY={binary_dir}\n"
                # Replace in command with $BINARY/name
                command_copy[binary_idx] = f"$BINARY/{binary_name}"
        
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
