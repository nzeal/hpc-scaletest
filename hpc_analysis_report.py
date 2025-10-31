#!/usr/bin/env python3
"""
HPC Analysis Report - Extract timing data from completed jobs and generate scaling reports.

This script iterates through job output directories, extracts timing information from
job.out files, and generates comprehensive scaling analysis reports.

Usage:
    python hpc_analysis_report.py output/iPIC3D-CPU-NS_weak_20251025_174137
    python hpc_analysis_report.py output/iPIC3D-CPU-NS_weak_20251025_174137 --scaling weak
    python hpc_analysis_report.py output/iPIC3D-CPU-NS_weak_20251025_174137 --verbose
"""

import argparse
import sys
import logging
import json
import re
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class JobAnalyzer:
    """Analyze individual job output files."""
    
    def __init__(self, job_dir: Path):
        self.job_dir = job_dir
        self.job_name = job_dir.name
        self.job_out = job_dir / "job.out"
        self.job_err = job_dir / "job.err"
        self.metadata = job_dir / "metadata.json"
    
    def analyze(self) -> Dict:
        """Analyze job outputs and extract metrics."""
        result = {
            'job_name': self.job_name,
            'status': 'unknown',
            'timing': None,
            'errors': [],
            'warnings': [],
            'metadata': {}
        }
        
        # Read metadata if available
        if self.metadata.exists():
            try:
                with open(self.metadata, 'r') as f:
                    result['metadata'] = json.load(f)
            except Exception as e:
                logger.warning(f"Could not read metadata for {self.job_name}: {e}")
        
        # Check if output files exist
        if not self.job_out.exists():
            result['status'] = 'no_output'
            result['errors'].append(f"job.out not found - job may not have run")
            return result
        
        # Parse job.out for timing and errors
        timing = self._extract_timing()
        errors = self._extract_errors()
        
        if timing:
            result['timing'] = timing
            result['status'] = 'completed' if not errors else 'completed_with_errors'
        elif errors:
            result['status'] = 'failed'
        else:
            result['status'] = 'no_timing'
            result['warnings'].append("No timing data found in output")
        
        result['errors'].extend(errors)
        
        return result
    
    def _extract_timing(self) -> Optional[float]:
        """Extract execution time from job.out."""
        try:
            content = self.job_out.read_text()
            
            # Look for common timing patterns
            patterns = [
                r'Total of ([\d.]+) seconds elapsed for process',
                r'Total execution time:\s*([\d.]+)\s*seconds',
                r'Wall time:\s*([\d.]+)\s*s',
                r'Elapsed time:\s*([\d.]+)\s*seconds',
                r'Time:\s*([\d.]+)\s*s'
            ]
            
            for pattern in patterns:
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    timing = float(match.group(1))
                    logger.debug(f"{self.job_name}: Found timing = {timing}s")
                    return timing
            
            logger.warning(f"{self.job_name}: No timing pattern matched")
            return None
            
        except Exception as e:
            logger.error(f"Error extracting timing from {self.job_name}: {e}")
            return None
    
    def _extract_errors(self) -> List[str]:
        """Extract errors from job.out and job.err."""
        errors = []
        
        # Check job.err
        if self.job_err.exists():
            try:
                err_content = self.job_err.read_text()
                if err_content.strip():
                    # Look for actual errors (not just warnings)
                    error_lines = []
                    for line in err_content.splitlines():
                        line_lower = line.lower()
                        if any(keyword in line_lower for keyword in ['error:', 'failed', 'fatal', 'segmentation fault', 'abort']):
                            error_lines.append(line.strip())
                    
                    if error_lines:
                        errors.extend(error_lines)
            except Exception as e:
                logger.warning(f"Could not read job.err for {self.job_name}: {e}")
        
        # Check job.out for errors
        if self.job_out.exists():
            try:
                out_content = self.job_out.read_text()
                # Common error patterns
                error_patterns = [
                    r'cd:.*No such file or directory',
                    r'srun: error:.*',
                    r'sbatch: error:.*',
                    r'.*: command not found',
                    r'Permission denied',
                    r'Cannot allocate memory',
                    r'Segmentation fault'
                ]
                
                for pattern in error_patterns:
                    matches = re.findall(pattern, out_content, re.IGNORECASE)
                    errors.extend(matches)
            except Exception as e:
                logger.warning(f"Could not scan job.out for errors in {self.job_name}: {e}")
        
        return errors


class ScalingAnalyzer:
    """Analyze scaling performance across multiple jobs."""
    
    def __init__(self, run_dir: Path, scaling_type: str = 'weak'):
        self.run_dir = run_dir
        self.scaling_type = scaling_type
        self.jobs = {}
    
    def analyze_all_jobs(self):
        """Find and analyze all job directories."""
        logger.info(f"Scanning for job directories in: {self.run_dir}")
        
        job_dirs = sorted([d for d in self.run_dir.iterdir() 
                          if d.is_dir() and d.name.startswith('nodes_')])
        
        if not job_dirs:
            logger.error("No job directories found (expected: nodes_1, nodes_2, etc.)")
            return
        
        logger.info(f"Found {len(job_dirs)} job directories")
        
        for job_dir in job_dirs:
            logger.info(f"Analyzing: {job_dir.name}")
            analyzer = JobAnalyzer(job_dir)
            result = analyzer.analyze()
            self.jobs[job_dir.name] = result
            
            # Log job status
            status_icon = {
                'completed': '✓',
                'completed_with_errors': '⚠',
                'failed': '✗',
                'no_output': '?',
                'no_timing': '⊘',
                'unknown': '?'
            }.get(result['status'], '?')
            
            status_msg = f"  {status_icon} {job_dir.name}: {result['status']}"
            if result['timing']:
                status_msg += f" ({result['timing']:.2f}s)"
            logger.info(status_msg)
            
            if result['errors']:
                for error in result['errors'][:3]:  # Show first 3 errors
                    logger.warning(f"    Error: {error[:100]}")
    
    def generate_report(self) -> Tuple[str, Dict]:
        """Generate scaling analysis report."""
        report_lines = []
        report_lines.append("=" * 70)
        report_lines.append(f"HPC SCALING ANALYSIS REPORT")
        report_lines.append("=" * 70)
        report_lines.append(f"Run Directory: {self.run_dir.name}")
        report_lines.append(f"Scaling Type: {self.scaling_type.upper()}")
        report_lines.append(f"Analysis Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append(f"Total Jobs: {len(self.jobs)}")
        report_lines.append("=" * 70)
        report_lines.append("")
        
        # Extract successful jobs with timing
        successful_jobs = {name: data for name, data in self.jobs.items() 
                          if data['timing'] is not None}
        
        if not successful_jobs:
            report_lines.append("⚠ NO TIMING DATA FOUND")
            report_lines.append("")
            report_lines.append("Possible reasons:")
            report_lines.append("  - Jobs are still running")
            report_lines.append("  - Jobs failed before completion")
            report_lines.append("  - Output files are missing or corrupted")
            report_lines.append("")
            
            # Show job statuses
            report_lines.append("Job Statuses:")
            for name, data in sorted(self.jobs.items()):
                report_lines.append(f"  {name}: {data['status']}")
                if data['errors']:
                    for error in data['errors'][:2]:
                        report_lines.append(f"    - {error[:80]}")
        else:
            # Timing results
            report_lines.append("TIMING RESULTS:")
            report_lines.append("-" * 70)
            for name in sorted(successful_jobs.keys(), key=lambda x: int(x.split('_')[1])):
                timing = successful_jobs[name]['timing']
                nodes = int(name.split('_')[1])
                report_lines.append(f"  {name:15s}: {timing:10.2f} seconds ({nodes:3d} nodes)")
            report_lines.append("")
            
            # Calculate scaling metrics
            if self.scaling_type == 'weak':
                metrics = self._calculate_weak_scaling(successful_jobs)
            else:
                metrics = self._calculate_strong_scaling(successful_jobs)
            
            report_lines.extend(metrics['lines'])
        
        report_lines.append("")
        report_lines.append("=" * 70)
        report_lines.append("END OF REPORT")
        report_lines.append("=" * 70)
        
        report_text = "\n".join(report_lines)
        
        # Generate JSON report
        json_report = {
            'run_directory': str(self.run_dir),
            'scaling_type': self.scaling_type,
            'analysis_time': datetime.now().isoformat(),
            'total_jobs': len(self.jobs),
            'successful_jobs': len(successful_jobs),
            'jobs': self.jobs
        }
        
        return report_text, json_report
    
    def _calculate_weak_scaling(self, jobs: Dict) -> Dict:
        """Calculate weak scaling efficiency."""
        lines = []
        lines.append("WEAK SCALING EFFICIENCY:")
        lines.append("-" * 70)
        
        # Sort by node count
        sorted_jobs = sorted(jobs.items(), key=lambda x: int(x[0].split('_')[1]))
        
        if not sorted_jobs:
            lines.append("  No data available")
            return {'lines': lines}
        
        # Use first job as baseline
        baseline_name, baseline_data = sorted_jobs[0]
        baseline_time = baseline_data['timing']
        baseline_nodes = int(baseline_name.split('_')[1])
        
        lines.append(f"  Baseline: {baseline_name} = {baseline_time:.2f}s (100.00%)")
        
        efficiencies = []
        for name, data in sorted_jobs[1:]:
            timing = data['timing']
            nodes = int(name.split('_')[1])
            efficiency = (baseline_time / timing) * 100
            efficiencies.append(efficiency)
            lines.append(f"  {name:15s}: {timing:8.2f}s ({efficiency:6.2f}% efficiency)")
        
        if efficiencies:
            avg_efficiency = sum(efficiencies) / len(efficiencies)
            lines.append("")
            lines.append(f"  Average Efficiency (excluding baseline): {avg_efficiency:.2f}%")
        
        lines.append("")
        return {'lines': lines, 'efficiencies': efficiencies}
    
    def _calculate_strong_scaling(self, jobs: Dict) -> Dict:
        """Calculate strong scaling speedup and efficiency."""
        lines = []
        lines.append("STRONG SCALING ANALYSIS:")
        lines.append("-" * 70)
        
        # Sort by node count
        sorted_jobs = sorted(jobs.items(), key=lambda x: int(x[0].split('_')[1]))
        
        if not sorted_jobs:
            lines.append("  No data available")
            return {'lines': lines}
        
        # Use first job as baseline
        baseline_name, baseline_data = sorted_jobs[0]
        baseline_time = baseline_data['timing']
        baseline_nodes = int(baseline_name.split('_')[1])
        
        lines.append(f"  Baseline: {baseline_name} = {baseline_time:.2f}s")
        lines.append("")
        lines.append(f"  {'Job':<15s} {'Time (s)':<12s} {'Speedup':<10s} {'Efficiency':<12s}")
        lines.append("  " + "-" * 60)
        
        speedups = []
        efficiencies = []
        
        for name, data in sorted_jobs:
            timing = data['timing']
            nodes = int(name.split('_')[1])
            speedup = baseline_time / timing
            efficiency = (speedup / (nodes / baseline_nodes)) * 100
            
            speedups.append(speedup)
            efficiencies.append(efficiency)
            
            lines.append(f"  {name:<15s} {timing:<12.2f} {speedup:<10.2f}x {efficiency:<12.2f}%")
        
        if len(speedups) > 1:
            avg_efficiency = sum(efficiencies[1:]) / len(efficiencies[1:])
            lines.append("")
            lines.append(f"  Average Parallel Efficiency: {avg_efficiency:.2f}%")
        
        lines.append("")
        return {'lines': lines, 'speedups': speedups, 'efficiencies': efficiencies}


def main():
    parser = argparse.ArgumentParser(
        description='Analyze HPC job outputs and generate scaling reports',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        'run_directory',
        type=Path,
        help='Path to the test run directory (e.g., output/iPIC3D-CPU-NS_weak_20251025_174137)'
    )
    
    parser.add_argument(
        '--scaling',
        choices=['weak', 'strong'],
        help='Scaling type (auto-detected from directory name if not specified)'
    )
    
    parser.add_argument(
        '--output',
        type=Path,
        help='Output file for the report (default: <run_dir>/AnalysisReport.txt)'
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
    
    # Auto-detect scaling type
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
    
    logger.info(f"{'='*70}")
    logger.info(f"HPC ANALYSIS REPORT")
    logger.info(f"{'='*70}")
    logger.info(f"Run directory: {run_dir}")
    logger.info(f"Scaling type: {scaling_type}")
    logger.info(f"{'='*70}\n")
    
    # Analyze jobs
    analyzer = ScalingAnalyzer(run_dir, scaling_type)
    analyzer.analyze_all_jobs()
    
    # Generate reports
    logger.info("\nGenerating reports...")
    report_text, json_report = analyzer.generate_report()
    
    # Save text report
    if args.output:
        report_path = args.output
    else:
        report_path = run_dir / "AnalysisReport.txt"
    
    with open(report_path, 'w') as f:
        f.write(report_text)
    
    # Save JSON report
    json_path = run_dir / f"{scaling_type}_analysis_report.json"
    with open(json_path, 'w') as f:
        json.dump(json_report, f, indent=2)
    
    logger.info(f"\n{'='*70}")
    logger.info(f"✅ ANALYSIS COMPLETE")
    logger.info(f"{'='*70}")
    logger.info(f"Text report: {report_path}")
    logger.info(f"JSON report: {json_path}")
    logger.info(f"{'='*70}\n")
    
    # Print report to console
    logger.info("Report Preview:")
    logger.info("-" * 70)
    print(report_text)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
