# backends/schedulers/slurm.py
import subprocess
import re
import os
import time
from pathlib import Path
from typing import Dict, Any, Optional
from core.types import JobID, JobStatus
from core.abstracts import AbstractScheduler
import logging

logger = logging.getLogger(__name__)

class SlurmScheduler(AbstractScheduler):
    _STATE_MAP_SQUEUE = {
        'pd': JobStatus.PENDING,
        'r': JobStatus.RUNNING,
        'cg': JobStatus.RUNNING,  # Completing, treat as running
        'cf': JobStatus.RUNNING,  # Configuring
        's': JobStatus.PENDING,   # Suspended
        'st': JobStatus.PENDING,  # Stopped
    }

    _STATE_MAP_SACCT = {
        'completed': JobStatus.COMPLETED,
        'cd': JobStatus.COMPLETED,
        'failed': JobStatus.FAILED,
        'f': JobStatus.FAILED,
        'cancelled': JobStatus.FAILED,
        'ca': JobStatus.FAILED,
        'timeout': JobStatus.FAILED,
        'to': JobStatus.FAILED,
        'out_of_memory': JobStatus.FAILED,
        'oom': JobStatus.FAILED,
        'node_fail': JobStatus.FAILED,
        'nf': JobStatus.FAILED,
        'boot_fail': JobStatus.FAILED,
        'bf': JobStatus.FAILED,
        'preempted': JobStatus.FAILED,
        'pr': JobStatus.FAILED,
        'deadline': JobStatus.FAILED,
        'dl': JobStatus.FAILED,
    }

    def submit(self, job_script: Path) -> JobID:
        """Submit a job script using sbatch and return the job ID."""
        try:
            result = subprocess.run(
                ['sbatch', str(job_script)],
                capture_output=True,
                text=True,
                check=True,
                env=os.environ
            )
            job_id_match = re.search(r'Submitted batch job (\d+)', result.stdout)
            if not job_id_match:
                raise RuntimeError(f"Could not extract job ID from sbatch output: {result.stdout}")
            job_id = JobID(job_id_match.group(1))
            logger.info(f"Submitted Slurm job {job_id}")
            return job_id
        except subprocess.CalledProcessError as e:
            logger.error(f"Slurm submission failed: {e.stderr}")
            raise RuntimeError(f"Slurm submission failed: {e.stderr}")

    def monitor(self, job_id: JobID) -> JobStatus:
        """Monitor the job status using squeue first, fallback to sacct."""
        # Try squeue for active jobs
        try:
            result = subprocess.run(
                ['squeue', '-j', str(job_id), '--noheader', '--format=%T'],
                capture_output=True,
                text=True,
                check=False,
                env=os.environ
            )
            status_str = result.stdout.strip().lower()
            if status_str:
                # Map squeue abbreviated state
                for abbr, mapped_status in self._STATE_MAP_SQUEUE.items():
                    if abbr in status_str:
                        return mapped_status
                # If not mapped, assume pending or running based on common
                if 'pending' in status_str or status_str == 'pd':
                    return JobStatus.PENDING
                elif 'running' in status_str or status_str == 'r':
                    return JobStatus.RUNNING
                else:
                    logger.warning(f"Unknown squeue state: {status_str}")
                    return JobStatus.RUNNING  # Default to running if active
        except Exception as e:
            logger.debug(f"squeue failed for job {job_id}: {e}")

        # Fallback to sacct for completed/failed jobs
        try:
            sacct_result = subprocess.run(
                ['sacct', '-j', str(job_id), '--format=State', '--noheader', '--parsable2'],
                capture_output=True,
                text=True,
                check=False,
                env=os.environ
            )
            state_str = sacct_result.stdout.strip().lower()
            if state_str:
                for full_state, mapped_status in self._STATE_MAP_SACCT.items():
                    if full_state in state_str:
                        return mapped_status
                logger.warning(f"Unknown sacct state: {state_str}")
                return JobStatus.COMPLETED  # Assume completed if not matched
        except Exception as e:
            logger.debug(f"sacct failed for job {job_id}: {e}")

        # If neither, assume completed
        return JobStatus.COMPLETED

    def cancel(self, job_id: JobID) -> bool:
        """Cancel the job using scancel."""
        try:
            result = subprocess.run(
                ['scancel', str(job_id)],
                capture_output=True,
                check=True,
                env=os.environ
            )
            logger.info(f"Cancelled Slurm job {job_id}")
            return True
        except subprocess.CalledProcessError:
            logger.error(f"Failed to cancel Slurm job {job_id}")
            return False

    def wait(self, job_id: JobID, timeout: Optional[int] = None) -> Dict[str, Any]:
        """Wait for job completion, polling status, and gather results using sacct and output file."""
        start_time = time.time()
        poll_interval = 30  # Poll every 30 seconds
        while True:
            status = self.monitor(job_id)
            if status in (JobStatus.COMPLETED, JobStatus.FAILED):
                break
            elapsed = time.time() - start_time
            if timeout and elapsed > timeout:
                logger.warning(f"Timeout waiting for Slurm job {job_id}")
                self.cancel(job_id)
                return {'status': JobStatus.FAILED, 'error': 'Timeout'}
            time.sleep(poll_interval)

        # Gather results
        results = {'status': status, 'job_id': str(job_id)}

        # Use sacct for accounting info
        try:
            sacct_result = subprocess.run(
                ['sacct', '-j', str(job_id), '--format=JobID,JobName,Partition,State,ExitCode,Elapsed,MaxRSS,CPUTime,NTasks,AllocCPUS,ReqMem,NodeList',
                 '--noheader', '--parsable2'],
                capture_output=True,
                text=True,
                check=False,
                env=os.environ
            )
            if sacct_result.stdout.strip():
                # Split by | for parsable2
                details = sacct_result.stdout.strip().split('|')
                results['sacct_details'] = dict(zip(['JobID', 'JobName', 'Partition', 'State', 'ExitCode', 'Elapsed', 'MaxRSS', 'CPUTime', 'NTasks', 'AllocCPUS', 'ReqMem', 'NodeList'], details))
            else:
                logger.warning(f"No sacct output for job {job_id}")
        except Exception as e:
            logger.error(f"sacct failed for job {job_id}: {e}")

        # Read output file for stdout/stderr and custom metrics
        out_file = Path(f"slurm-{job_id}.out")
        err_file = Path(f"slurm-{job_id}.err")
        if out_file.exists():
            results['stdout'] = out_file.read_text()
        if err_file.exists():
            results['stderr'] = err_file.read_text()

        # Optionally parse runtime from stdout if benchmark outputs it
        # Assume benchmark outputs "Runtime: Xs" or similar
        if 'stdout' in results:
            runtime_match = re.search(r'Runtime:\s*(\d+\.?\d*)s?', results['stdout'], re.IGNORECASE)
            if runtime_match:
                results['parsed_runtime'] = float(runtime_match.group(1))

        logger.info(f"Slurm job {job_id} completed with status {status}")
        return results
