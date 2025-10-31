#!/usr/bin/env python3
"""
Utility functions for submitting prepared job scripts.
"""

import json
import subprocess
import re
import logging
from pathlib import Path
from typing import List, Tuple, Optional

logger = logging.getLogger(__name__)


def submit_prepared_jobs(output_dir: Path) -> Tuple[List[Tuple[str, str]], List[str]]:
    """
    Submit all prepared job scripts in the output directory.
    
    Args:
        output_dir: Path to the output directory containing job folders
        
    Returns:
        Tuple of (successful_submissions, failed_submissions) where:
        - successful_submissions: List of (job_id, slurm_job_id) tuples
        - failed_submissions: List of job_ids that failed to submit
    """
    successful_submissions = []
    failed_submissions = []
    
    # Find all metadata files
    metadata_files = list(output_dir.rglob("metadata.json"))
    
    for metadata_file in metadata_files:
        job_id = 'unknown'
        try:
            # Read metadata
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
            
            # Check if this job was prepared but not submitted
            if metadata.get('status') == 'prepared':
                job_id = metadata.get('job_id', 'unknown')
                script_path = Path(metadata.get('script_path', ''))
                
                if not script_path.exists():
                    logger.error(f"Job script not found for {job_id}: {script_path}")
                    failed_submissions.append(job_id)
                    continue
                
                # Submit the job
                logger.info(f"Submitting prepared job {job_id} using script {script_path}")
                result = subprocess.run(
                    ["sbatch", str(script_path)],
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                # Parse job ID
                match = re.search(r'Submitted batch job (\d+)', result.stdout)
                if match:
                    slurm_job_id = match.group(1)
                    logger.info(f"Successfully submitted job {job_id} as Slurm job {slurm_job_id}")
                    
                    # Update metadata
                    metadata['status'] = 'submitted'
                    metadata['submitted_at'] = __import__('datetime').datetime.now().isoformat()
                    metadata['slurm_job_id'] = slurm_job_id
                    
                    with open(metadata_file, 'w') as f:
                        json.dump(metadata, f, indent=2)
                    
                    successful_submissions.append((job_id, slurm_job_id))
                else:
                    logger.error(f"Could not parse job ID from sbatch output for {job_id}: {result.stdout}")
                    failed_submissions.append(job_id)
                    
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to submit job {job_id}: {e.stderr}")
            failed_submissions.append(job_id)
        except Exception as e:
            logger.error(f"Unexpected error submitting job {job_id}: {e}")
            failed_submissions.append(job_id)
    
    return successful_submissions, failed_submissions


def submit_single_job(script_path: Path) -> Optional[str]:
    """
    Submit a single job script.
    
    Args:
        script_path: Path to the job script (.sh file)
        
    Returns:
        Slurm job ID if successful, None otherwise
    """
    if not script_path.exists():
        logger.error(f"Job script not found: {script_path}")
        return None
    
    try:
        logger.info(f"Submitting job script: {script_path}")
        result = subprocess.run(
            ["sbatch", str(script_path)],
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse job ID
        match = re.search(r'Submitted batch job (\d+)', result.stdout)
        if match:
            job_id = match.group(1)
            logger.info(f"Successfully submitted job: {job_id}")
            return job_id
        else:
            logger.error(f"Could not parse job ID from sbatch output: {result.stdout}")
            return None
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to submit job: {e.stderr}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error submitting job: {e}")
        return None


def check_sbatch_availability() -> bool:
    """
    Check if sbatch is available in the current environment.
    
    Returns:
        True if sbatch is available, False otherwise
    """
    try:
        result = subprocess.run(["which", "sbatch"], capture_output=True, text=True)
        return result.returncode == 0
    except Exception:
        return False