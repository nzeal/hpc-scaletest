"""
Efficiency report generator for HPC scaling tests.
"""

import logging
import json
from pathlib import Path
from typing import List, Dict, Optional
from dataclasses import dataclass
from datetime import datetime

logger = logging.getLogger(__name__)


@dataclass
class ScalingResult:
    """Single scaling test result."""
    nodes: int
    procs: int
    time_seconds: float
    speedup: float = 0.0
    efficiency: float = 0.0
    procs_decomposition: tuple = (0, 0, 0)


class ReportGenerator:
    """Generates scaling efficiency reports from test results."""
    
    def __init__(self, output_dir: Path):
        """
        Initialize report generator.
        
        Args:
            output_dir: Directory containing test results
        """
        self.output_dir = output_dir
    
    def generate_scaling_report(
        self,
        scaling_type: str = "strong",
        output_file: Optional[str] = None
    ) -> str:
        """
        Generate a scaling efficiency report.
        
        Args:
            scaling_type: Type of scaling test ("strong" or "weak")
            output_file: Optional output filename (auto-generated if None)
            
        Returns:
            Path to generated report file
        """
        # Read summary.json
        summary_path = self.output_dir / "summary.json"
        if not summary_path.exists():
            raise FileNotFoundError(f"Summary file not found: {summary_path}")
        
        with open(summary_path, 'r') as f:
            summary = json.load(f)
        
        # Parse results
        results = self._parse_results(summary)
        
        # Calculate metrics
        if scaling_type.lower() == "strong":
            results = self._calculate_strong_scaling_metrics(results)
        else:
            results = self._calculate_weak_scaling_metrics(results)
        
        # Generate report
        report_content = self._format_report(results, scaling_type, summary)
        
        # Save report
        if output_file is None:
            output_file = f"{scaling_type.capitalize()}ScalingReport.txt"
        report_path = self.output_dir / output_file
        
        with open(report_path, 'w') as f:
            f.write(report_content)
        
        logger.info(f"Generated {scaling_type} scaling report: {report_path}")
        
        # Also generate JSON version for programmatic access
        json_report_path = self.output_dir / f"{scaling_type}_scaling_report.json"
        self._save_json_report(results, scaling_type, summary, json_report_path)
        
        return str(report_path)
    
    def _parse_results(self, summary: Dict) -> List[ScalingResult]:
        """Parse results from summary JSON."""
        results = []
        
        for job in summary.get('jobs', []):
            # Extract time from metrics
            time_seconds = None
            if 'wall_time' in job:
                time_seconds = job['wall_time']
            elif 'metrics' in job and 'wall_time' in job['metrics']:
                time_seconds = job['metrics']['wall_time']
            
            if time_seconds is None:
                logger.warning(f"No timing data for job {job['job_id']}, skipping")
                continue
            
            result = ScalingResult(
                nodes=job['num_nodes'],
                procs=job['num_procs'],
                time_seconds=float(time_seconds),
                procs_decomposition=tuple(job.get('procs_decomposition', (0, 0, 0)))
            )
            results.append(result)
        
        # Sort by number of nodes
        results.sort(key=lambda x: x.nodes)
        
        return results
    
    def _calculate_strong_scaling_metrics(
        self,
        results: List[ScalingResult]
    ) -> List[ScalingResult]:
        """Calculate speedup and efficiency for strong scaling."""
        if not results:
            return results
        
        # Use first result as baseline
        baseline_time = results[0].time_seconds
        baseline_procs = results[0].procs
        
        for result in results:
            # Speedup = T_baseline / T_current
            result.speedup = baseline_time / result.time_seconds
            
            # Efficiency = Speedup / (Procs_current / Procs_baseline)
            proc_ratio = result.procs / baseline_procs
            result.efficiency = (result.speedup / proc_ratio) * 100.0
        
        return results
    
    def _calculate_weak_scaling_metrics(
        self,
        results: List[ScalingResult]
    ) -> List[ScalingResult]:
        """Calculate efficiency for weak scaling."""
        if not results:
            return results
        
        # Use first result as baseline
        baseline_time = results[0].time_seconds
        
        for result in results:
            # For weak scaling, ideal time should remain constant
            # Efficiency = T_baseline / T_current * 100%
            result.efficiency = (baseline_time / result.time_seconds) * 100.0
            
            # Speedup is not as meaningful for weak scaling, but we can show it
            result.speedup = baseline_time / result.time_seconds
        
        return results
    
    def _format_report(
        self,
        results: List[ScalingResult],
        scaling_type: str,
        summary: Dict
    ) -> str:
        """Format results as a text report."""
        lines = []
        
        # Header
        lines.append("=" * 80)
        lines.append(f"{scaling_type.capitalize()} Scaling Efficiency Report")
        lines.append("=" * 80)
        lines.append(f"Test Name: {summary.get('test_name', 'Unknown')}")
        lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        if 'completed_at' in summary:
            lines.append(f"Test Completed: {summary['completed_at']}")
        lines.append("=" * 80)
        lines.append("")
        
        # Table header
        lines.append(f"{'Nodes':<10} {'Procs':<10} {'Time(s)':<12} {'Speedup':<12} {'Efficiency':<12} {'Decomposition':<20}")
        lines.append("-" * 80)
        
        # Data rows
        for result in results:
            decomp_str = f"({result.procs_decomposition[0]}×{result.procs_decomposition[1]}×{result.procs_decomposition[2]})"
            lines.append(
                f"{result.nodes:<10} "
                f"{result.procs:<10} "
                f"{result.time_seconds:<12.2f} "
                f"{result.speedup:<12.2f} "
                f"{result.efficiency:<12.1f}% "
                f"{decomp_str:<20}"
            )
        
        lines.append("=" * 80)
        lines.append("")
        
        # Summary statistics
        lines.append("Summary Statistics:")
        lines.append("-" * 40)
        if results:
            avg_efficiency = sum(r.efficiency for r in results) / len(results)
            max_efficiency = max(r.efficiency for r in results)
            min_efficiency = min(r.efficiency for r in results)
            
            lines.append(f"  Average Efficiency: {avg_efficiency:.1f}%")
            lines.append(f"  Maximum Efficiency: {max_efficiency:.1f}%")
            lines.append(f"  Minimum Efficiency: {min_efficiency:.1f}%")
            lines.append("")
            
            if scaling_type.lower() == "strong":
                max_speedup = max(r.speedup for r in results)
                parallel_efficiency = (max_speedup / results[-1].procs) * results[0].procs * 100
                lines.append(f"  Maximum Speedup: {max_speedup:.2f}x")
                lines.append(f"  Parallel Efficiency (max config): {parallel_efficiency:.1f}%")
        
        lines.append("=" * 80)
        
        # Interpretation
        lines.append("")
        lines.append("Interpretation Guide:")
        lines.append("-" * 40)
        if scaling_type.lower() == "strong":
            lines.append("  Strong Scaling: Problem size is fixed, parallelism increases")
            lines.append("  - Speedup: How much faster compared to baseline configuration")
            lines.append("  - Efficiency: How well resources are utilized (100% = ideal)")
            lines.append("  - Good strong scaling: Efficiency > 70% up to high node counts")
        else:
            lines.append("  Weak Scaling: Problem size grows with parallelism")
            lines.append("  - Efficiency: How well execution time stays constant (100% = ideal)")
            lines.append("  - Good weak scaling: Efficiency > 80% across all configurations")
        lines.append("=" * 80)
        
        return "\n".join(lines)
    
    def _save_json_report(
        self,
        results: List[ScalingResult],
        scaling_type: str,
        summary: Dict,
        output_path: Path
    ):
        """Save report in JSON format."""
        report_data = {
            'test_name': summary.get('test_name', 'Unknown'),
            'scaling_type': scaling_type,
            'generated_at': datetime.now().isoformat(),
            'completed_at': summary.get('completed_at'),
            'results': [
                {
                    'nodes': r.nodes,
                    'procs': r.procs,
                    'time_seconds': r.time_seconds,
                    'speedup': r.speedup,
                    'efficiency_percent': r.efficiency,
                    'procs_decomposition': list(r.procs_decomposition)
                }
                for r in results
            ]
        }
        
        # Add summary statistics
        if results:
            report_data['statistics'] = {
                'average_efficiency': sum(r.efficiency for r in results) / len(results),
                'max_efficiency': max(r.efficiency for r in results),
                'min_efficiency': min(r.efficiency for r in results),
            }
            
            if scaling_type.lower() == "strong":
                max_speedup = max(r.speedup for r in results)
                report_data['statistics']['max_speedup'] = max_speedup
        
        with open(output_path, 'w') as f:
            json.dump(report_data, f, indent=2)
        
        logger.info(f"Saved JSON report: {output_path}")
    
    def generate_comparison_report(
        self,
        test_dirs: List[Path],
        output_file: str = "comparison_report.txt"
    ) -> str:
        """
        Generate a comparison report across multiple test runs.
        
        Args:
            test_dirs: List of test output directories
            output_file: Output filename
            
        Returns:
            Path to generated report
        """
        all_results = {}
        
        for test_dir in test_dirs:
            summary_path = test_dir / "summary.json"
            if not summary_path.exists():
                logger.warning(f"Summary not found in {test_dir}, skipping")
                continue
            
            with open(summary_path, 'r') as f:
                summary = json.load(f)
            
            test_name = summary.get('test_name', test_dir.name)
            all_results[test_name] = self._parse_results(summary)
        
        # Generate comparison report
        lines = []
        lines.append("=" * 100)
        lines.append("Multi-Test Comparison Report")
        lines.append("=" * 100)
        lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("")
        
        # Compare at each node count
        for test_name, results in all_results.items():
            lines.append(f"\n{test_name}:")
            lines.append("-" * 50)
            for result in results:
                lines.append(
                    f"  {result.nodes} nodes: {result.time_seconds:.2f}s "
                    f"(Efficiency: {result.efficiency:.1f}%)"
                )
        
        lines.append("\n" + "=" * 100)
        
        report_content = "\n".join(lines)
        output_path = self.output_dir / output_file
        
        with open(output_path, 'w') as f:
            f.write(report_content)
        
        logger.info(f"Generated comparison report: {output_path}")
        return str(output_path)
