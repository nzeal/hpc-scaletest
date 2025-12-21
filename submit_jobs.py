#!/usr/bin/env python3
"""
Job Submission Utility for HPC-ScaleTest
Submits prepared job scripts from a test output directory.
"""

import argparse
import sys
from pathlib import Path
from typing import List, Dict
import logging

from utils.job_submitter import submit_single_job
from utils.logging_config import setup_logging


def find_job_scripts(output_dir: Path) -> List[Path]:
    """
    Find all job scripts in the output directory.
    
    Args:
        output_dir: Output directory containing job scripts
        
    Returns:
        List of job script paths
    """
    job_scripts = []
    
    # Look for job.sh files in node subdirectories
    for node_dir in sorted(output_dir.glob("nodes_*")):
        if node_dir.is_dir():
            job_script = node_dir / "job.sh"
            if job_script.exists():
                job_scripts.append(job_script)
    
    # Also check for job scripts directly in output directory
    for job_script in output_dir.glob("*.sh"):
        if job_script not in job_scripts:
            job_scripts.append(job_script)
    
    return sorted(job_scripts)


def submit_jobs(output_dir: Path, dry_run: bool = False) -> Dict[str, str]:
    """
    Submit all job scripts in the output directory.
    
    Args:
        output_dir: Output directory containing job scripts
        dry_run: If True, only show what would be submitted
        
    Returns:
        Dictionary mapping job script paths to job IDs
    """
    logger = logging.getLogger(__name__)
    
    job_scripts = find_job_scripts(output_dir)
    
    if not job_scripts:
        logger.error(f"No job scripts found in {output_dir}")
        return {}
    
    logger.info(f"Found {len(job_scripts)} job script(s) to submit")
    
    submitted_jobs = {}
    failed_jobs = []
    
    for job_script in job_scripts:
        try:
            if dry_run:
                logger.info(f"[DRY RUN] Would submit: {job_script}")
                submitted_jobs[str(job_script)] = "DRY_RUN"
            else:
                logger.info(f"Submitting: {job_script}")
                job_id = submit_single_job(job_script)
                submitted_jobs[str(job_script)] = job_id
                logger.info(f"  → Job ID: {job_id}")
        except Exception as e:
            logger.error(f"Failed to submit {job_script}: {e}")
            failed_jobs.append(job_script)
    
    # Summary
    logger.info(f"\n{'='*60}")
    logger.info(f"Submission Summary:")
    logger.info(f"  Total scripts: {len(job_scripts)}")
    logger.info(f"  Successfully submitted: {len(submitted_jobs)}")
    logger.info(f"  Failed: {len(failed_jobs)}")
    
    if failed_jobs:
        logger.error(f"\nFailed job scripts:")
        for job_script in failed_jobs:
            logger.error(f"  - {job_script}")
    
    return submitted_jobs


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Submit prepared job scripts from HPC-ScaleTest output directory",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Submit all jobs in output directory
  python submit_jobs.py output/my_test_strong_20241209_120000

  # Dry run to see what would be submitted
  python submit_jobs.py output/my_test_strong_20241209_120000 --dry-run

  # Submit with verbose logging
  python submit_jobs.py output/my_test_strong_20241209_120000 --verbose
        """
    )
    
    parser.add_argument(
        'output_dir',
        type=Path,
        help='Output directory containing job scripts'
    )
    
    parser.add_argument(
        '--dry-run', '-d',
        action='store_true',
        help='Show what would be submitted without actually submitting'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    setup_logging(log_level)
    logger = logging.getLogger(__name__)
    
    # Validate output directory
    if not args.output_dir.exists():
        logger.error(f"Output directory does not exist: {args.output_dir}")
        sys.exit(1)
    
    if not args.output_dir.is_dir():
        logger.error(f"Path is not a directory: {args.output_dir}")
        sys.exit(1)
    
    # Submit jobs
    try:
        submitted_jobs = submit_jobs(args.output_dir, dry_run=args.dry_run)
        
        if not submitted_jobs:
            sys.exit(1)
        
        if not args.dry_run:
            # Write job IDs to file
            job_ids_file = args.output_dir / "job_ids.txt"
            with open(job_ids_file, 'w') as f:
                for job_script, job_id in submitted_jobs.items():
                    f.write(f"{job_id}\t{job_script}\n")
            logger.info(f"\nJob IDs written to: {job_ids_file}")
    
    except KeyboardInterrupt:
        logger.info("\nSubmission interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Submission failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
