#!/usr/bin/env python3
"""
Standalone report generator for HPC scaling tests.
Run this after jobs complete to generate performance reports.

Usage:
    python generate_report.py output/iPIC3D-CPU-NS_weak_20251025_163227
    python generate_report.py output/iPIC3D-CPU-NS_weak_20251025_163227 --scaling weak
"""

import argparse
import sys
import logging
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from utils.report_generator import ReportGenerator

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description='Generate HPC scaling reports from completed test runs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect scaling type from directory name
  python generate_report.py output/iPIC3D-CPU-NS_weak_20251025_163227
  
  # Explicitly specify scaling type
  python generate_report.py output/test_strong_20251025_120000 --scaling strong
  
  # Generate report with verbose output
  python generate_report.py output/test_weak_20251025_120000 --verbose
"""
    )
    
    parser.add_argument(
        'run_directory',
        type=Path,
        help='Path to the test run directory (e.g., output/iPIC3D-CPU-NS_weak_20251025_163227)'
    )
    
    parser.add_argument(
        '--scaling',
        choices=['weak', 'strong'],
        help='Scaling type (auto-detected from directory name if not specified)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validate run directory
    run_dir = args.run_directory.resolve()
    if not run_dir.exists():
        logger.error(f"Run directory not found: {run_dir}")
        sys.exit(1)
    
    if not run_dir.is_dir():
        logger.error(f"Not a directory: {run_dir}")
        sys.exit(1)
    
    # Auto-detect scaling type from directory name
    scaling_type = args.scaling
    if not scaling_type:
        dir_name = run_dir.name.lower()
        if '_weak_' in dir_name:
            scaling_type = 'weak'
        elif '_strong_' in dir_name:
            scaling_type = 'strong'
        else:
            logger.error("Could not auto-detect scaling type from directory name.")
            logger.error("Please specify --scaling weak or --scaling strong")
            sys.exit(1)
        logger.info(f"Auto-detected scaling type: {scaling_type}")
    
    logger.info(f"{'='*60}")
    logger.info(f"REPORT GENERATION")
    logger.info(f"{'='*60}")
    logger.info(f"Run directory: {run_dir}")
    logger.info(f"Scaling type: {scaling_type}")
    logger.info(f"{'='*60}\n")
    
    # Check for job output files
    job_dirs = [d for d in run_dir.iterdir() if d.is_dir() and d.name.startswith('nodes_')]
    if not job_dirs:
        logger.warning("No job directories found (expected: nodes_1, nodes_2, etc.)")
    
    logger.info(f"Found {len(job_dirs)} job directories")
    
    # Check for timing data
    timing_found = False
    for job_dir in job_dirs:
        out_file = job_dir / "job.out"
        if out_file.exists():
            timing_found = True
            logger.debug(f"  {job_dir.name}: job.out found")
        else:
            logger.warning(f"  {job_dir.name}: job.out not found (job may still be running)")
    
    if not timing_found:
        logger.warning("\n⚠ No job output files found!")
        logger.warning("Jobs may still be running or failed to produce output.")
        logger.info("\nTo check job status:")
        logger.info("  squeue -u $USER")
        logger.info("  sacct --format=JobID,JobName,State,Start,End")
        logger.info("\nTo check job output:")
        for job_dir in job_dirs[:2]:  # Show first 2 as examples
            logger.info(f"  ls -la {job_dir}/")
        logger.info("")
    
    # Generate report
    try:
        logger.info("Generating report...\n")
        generator = ReportGenerator(run_dir)
        report_path = generator.generate_scaling_report(scaling_type=scaling_type)
        
        logger.info(f"\n{'='*60}")
        logger.info("✅ REPORT GENERATED SUCCESSFULLY")
        logger.info(f"{'='*60}")
        logger.info(f"Text report: {report_path}")
        
        # Check if JSON report exists
        json_report = run_dir / f"{scaling_type}_scaling_report.json"
        if json_report.exists():
            logger.info(f"JSON report: {json_report}")
        
        logger.info(f"{'='*60}\n")
        
        # Show preview of report
        if report_path.exists():
            logger.info("Report preview (first 30 lines):")
            logger.info("-" * 60)
            with open(report_path, 'r') as f:
                lines = f.readlines()
                for line in lines[:30]:
                    print(line.rstrip())
                if len(lines) > 30:
                    logger.info(f"... ({len(lines) - 30} more lines)")
            logger.info("-" * 60)
        
        return 0
        
    except Exception as e:
        logger.error(f"\n❌ Report generation failed: {e}", exc_info=args.verbose)
        return 1


if __name__ == '__main__':
    sys.exit(main())
