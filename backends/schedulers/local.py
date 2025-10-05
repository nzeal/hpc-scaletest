# backends/schedulers/local.py
import subprocess
import os
from pathlib import Path
from typing import Dict, Any, Optional
from ...core.types import JobID, JobStatus
from ...core.abstracts import AbstractScheduler
import logging

logger = logging.getLogger(__name__)

class LocalScheduler(AbstractScheduler):
    def __init__(self):
        self._job_counter = 0

    def submit(self, job_script: Path) -> JobID:
        self._job_counter += 1
        job_id = f"local_{self._job_counter}"
        logger.info(f"Submitted local job {job_id} (script: {job_script})")
        # For local, execution happens in wait; just return ID
        return JobID(job_id)

    def monitor(self, job_id: JobID) -> JobStatus:
        # For local, assume immediately completed after submit
        logger.debug(f"Monitoring local job {job_id}: COMPLETED")
        return JobStatus.COMPLETED

    def cancel(self, job_id: JobID) -> bool:
        # No-op for local, but log
        logger.info(f"Cancel requested for local job {job_id} (no-op)")
        return True

    def wait(self, job_id: JobID, timeout: Optional[int] = None) -> Dict[str, Any]:
        # Execute the script synchronously with timeout
        try:
            result = subprocess.run(
                ['bash', str(job_script)],
                capture_output=True,
                text=True,
                timeout=timeout,
                env=os.environ
            )
            status = JobStatus.COMPLETED if result.returncode == 0 else JobStatus.FAILED
            logger.info(f"Local job {job_id} completed with status {status}")
            return {
                'status': status,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'returncode': result.returncode,
                'job_id': str(job_id)
            }
        except subprocess.TimeoutExpired:
            logger.warning(f"Local job {job_id} timed out")
            return {'status': JobStatus.FAILED, 'error': 'Timeout'}
        except Exception as e:
            logger.error(f"Local job {job_id} failed: {e}")
            return {'status': JobStatus.FAILED, 'error': str(e)}
