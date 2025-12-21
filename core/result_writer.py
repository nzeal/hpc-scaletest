#!/usr/bin/env python3
"""
core/result_writer.py - Incremental Result Writer

Writes scaling results incrementally as each run completes, rather than
waiting for all runs to finish.

Output Structure:
    data/
    ├── weak_scaling_report_node1.json
    ├── weak_scaling_report_node2.json
    ├── weak_scaling_report_node4.json
    ├── summary.json
    └── weak_scaling_report.json (aggregated)

Design:
- Each run's result is written immediately upon completion
- Results can be analyzed while runs are still in progress
- Final aggregation happens after all runs complete (optional)

Author: HPC-ScaleTest Contributors
"""

import json
import logging
import os
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, field, asdict
from typing import Optional, List, Dict, Any

logger = logging.getLogger(__name__)


@dataclass
class RunResult:
    """Result from a single scaling run."""
    
    # Run identification
    run_id: str
    scaling_type: str  # "weak" or "strong"
    node_count: int
    
    # Timing
    start_time: str = ""
    end_time: str = ""
    wall_time_seconds: float = 0.0
    
    # Resource configuration
    total_mpi_ranks: int = 0
    mpi_ranks_per_node: int = 0
    cores_per_rank: int = 0
    gpus_per_node: int = 0
    total_gpus: int = 0
    
    # MPI decomposition
    procs_x: int = 1
    procs_y: int = 1
    procs_z: int = 1
    
    # Problem configuration (for weak scaling)
    domain_x: float = 0.0
    domain_y: float = 0.0
    domain_z: float = 0.0
    cells_x: int = 0
    cells_y: int = 0
    cells_z: int = 0
    
    # Execution status
    status: str = "pending"  # pending, running, completed, failed
    exit_code: int = -1
    error_message: str = ""
    
    # Metrics (computed after completion)
    speedup: float = 0.0
    efficiency: float = 0.0
    
    # Paths
    job_directory: str = ""
    job_script: str = ""
    output_file: str = ""
    error_file: str = ""
    
    # SLURM info
    slurm_job_id: str = ""
    partition: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'RunResult':
        """Create from dictionary."""
        return cls(**data)


@dataclass 
class ScalingSummary:
    """Summary of all scaling runs."""
    
    # Test identification
    test_name: str
    scaling_type: str
    
    # Timing
    start_time: str = ""
    end_time: str = ""
    total_duration_seconds: float = 0.0
    
    # Configuration
    partition: str = ""
    cpus_per_node: int = 0
    gpus_per_node: int = 0
    
    # Run counts
    total_runs: int = 0
    completed_runs: int = 0
    failed_runs: int = 0
    
    # Node sequence
    node_sequence: List[int] = field(default_factory=list)
    
    # Baseline metrics
    baseline_time: float = 0.0
    baseline_nodes: int = 1
    
    # Results list
    results: List[RunResult] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        d = asdict(self)
        d['results'] = [r.to_dict() if isinstance(r, RunResult) else r for r in self.results]
        return d


class IncrementalResultWriter:
    """
    Writes scaling results incrementally as each run completes.
    
    Usage:
        writer = IncrementalResultWriter(output_dir, "weak", "boost_usr_prod")
        
        # After each run completes:
        writer.write_run_result(result)
        
        # After all runs complete:
        writer.write_summary()
        writer.write_aggregated_report()
    """
    
    def __init__(
        self,
        output_dir: Path,
        scaling_type: str,
        partition: str = "",
        test_name: str = "scaling_test"
    ):
        """
        Initialize the result writer.
        
        Args:
            output_dir: Base output directory for the test
            scaling_type: "weak" or "strong"
            partition: SLURM partition name
            test_name: Name of the test
        """
        self.output_dir = Path(output_dir)
        self.scaling_type = scaling_type.lower()
        self.partition = partition
        self.test_name = test_name
        
        # Create data directory
        self.data_dir = self.output_dir / "data"
        self.data_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize summary
        self.summary = ScalingSummary(
            test_name=test_name,
            scaling_type=scaling_type,
            partition=partition,
            start_time=datetime.now().isoformat()
        )
        
        # Track results
        self._results: Dict[int, RunResult] = {}  # node_count -> result
        self._baseline_time: Optional[float] = None
        
        logger.info(f"Initialized incremental result writer")
        logger.info(f"  Output directory: {self.data_dir}")
        logger.info(f"  Scaling type: {self.scaling_type}")
    
    def write_run_result(self, result: RunResult) -> Path:
        """
        Write a single run's result immediately.
        
        Args:
            result: RunResult from completed run
            
        Returns:
            Path to the written JSON file
        """
        # Generate filename
        filename = f"{self.scaling_type}_scaling_report_node{result.node_count}.json"
        filepath = self.data_dir / filename
        
        # Compute metrics if this run completed successfully
        if result.status == "completed" and result.wall_time_seconds > 0:
            self._compute_metrics(result)
        
        # Write JSON
        with open(filepath, 'w') as f:
            json.dump(result.to_dict(), f, indent=2)
        
        logger.info(f"✓ Wrote result: {filepath}")
        
        # Store in memory for summary
        self._results[result.node_count] = result
        
        # Update baseline if this is node 1
        if result.node_count == 1 and result.status == "completed":
            self._baseline_time = result.wall_time_seconds
            logger.info(f"  Baseline time set: {self._baseline_time:.2f}s")
        
        return filepath
    
    def _compute_metrics(self, result: RunResult) -> None:
        """Compute speedup and efficiency for a result."""
        if self._baseline_time is None:
            # No baseline yet - can't compute
            return
        
        if self._baseline_time <= 0:
            return
        
        if self.scaling_type == "strong":
            # Strong scaling: speedup = T1 / Tn
            result.speedup = self._baseline_time / result.wall_time_seconds
            # Efficiency = speedup / N
            result.efficiency = result.speedup / result.node_count * 100
        else:
            # Weak scaling: efficiency = T1 / Tn (ideally 100%)
            result.efficiency = (self._baseline_time / result.wall_time_seconds) * 100
            result.speedup = result.node_count  # Ideal speedup for weak scaling
    
    def update_run_status(
        self,
        node_count: int,
        status: str,
        wall_time: float = 0.0,
        exit_code: int = 0,
        error_message: str = ""
    ) -> None:
        """
        Update the status of an existing run and rewrite its JSON.
        
        Args:
            node_count: Number of nodes for this run
            status: New status ("running", "completed", "failed")
            wall_time: Wall clock time in seconds
            exit_code: Process exit code
            error_message: Error message if failed
        """
        if node_count not in self._results:
            logger.warning(f"No result found for node count {node_count}")
            return
        
        result = self._results[node_count]
        result.status = status
        result.wall_time_seconds = wall_time
        result.exit_code = exit_code
        result.error_message = error_message
        
        if status == "completed":
            result.end_time = datetime.now().isoformat()
        
        # Rewrite the JSON file
        self.write_run_result(result)
    
    def create_pending_result(
        self,
        node_count: int,
        total_mpi_ranks: int,
        mpi_ranks_per_node: int,
        cores_per_rank: int,
        gpus_per_node: int = 0,
        procs_decomp: tuple = (1, 1, 1),
        job_directory: str = "",
        slurm_job_id: str = ""
    ) -> RunResult:
        """
        Create a pending result before the run starts.
        
        Args:
            node_count: Number of nodes
            total_mpi_ranks: Total MPI processes
            mpi_ranks_per_node: MPI processes per node
            cores_per_rank: CPU cores per MPI process
            gpus_per_node: GPUs per node
            procs_decomp: MPI decomposition (px, py, pz)
            job_directory: Path to job directory
            slurm_job_id: SLURM job ID
            
        Returns:
            RunResult with pending status
        """
        result = RunResult(
            run_id=f"{self.scaling_type}_node{node_count}",
            scaling_type=self.scaling_type,
            node_count=node_count,
            start_time=datetime.now().isoformat(),
            total_mpi_ranks=total_mpi_ranks,
            mpi_ranks_per_node=mpi_ranks_per_node,
            cores_per_rank=cores_per_rank,
            gpus_per_node=gpus_per_node,
            total_gpus=node_count * gpus_per_node,
            procs_x=procs_decomp[0],
            procs_y=procs_decomp[1],
            procs_z=procs_decomp[2],
            status="pending",
            job_directory=job_directory,
            slurm_job_id=slurm_job_id,
            partition=self.partition
        )
        
        # Write immediately
        self.write_run_result(result)
        
        return result
    
    def write_summary(self) -> Path:
        """
        Write the summary.json file after all runs complete.
        
        Returns:
            Path to summary.json
        """
        self.summary.end_time = datetime.now().isoformat()
        
        # Calculate duration
        try:
            start = datetime.fromisoformat(self.summary.start_time)
            end = datetime.fromisoformat(self.summary.end_time)
            self.summary.total_duration_seconds = (end - start).total_seconds()
        except:
            pass
        
        # Count results
        self.summary.total_runs = len(self._results)
        self.summary.completed_runs = sum(1 for r in self._results.values() if r.status == "completed")
        self.summary.failed_runs = sum(1 for r in self._results.values() if r.status == "failed")
        
        # Extract node sequence
        self.summary.node_sequence = sorted(self._results.keys())
        
        # Set baseline
        if 1 in self._results:
            self.summary.baseline_time = self._results[1].wall_time_seconds
            self.summary.baseline_nodes = 1
        
        # Add results
        self.summary.results = list(self._results.values())
        
        # Write JSON
        filepath = self.data_dir / "summary.json"
        with open(filepath, 'w') as f:
            json.dump(self.summary.to_dict(), f, indent=2)
        
        logger.info(f"✓ Wrote summary: {filepath}")
        
        return filepath
    
    def write_aggregated_report(self) -> Path:
        """
        Write the aggregated scaling report after all runs complete.
        
        Returns:
            Path to aggregated report
        """
        # Sort results by node count
        sorted_results = sorted(self._results.values(), key=lambda r: r.node_count)
        
        # Build aggregated report
        report = {
            "test_name": self.test_name,
            "scaling_type": self.scaling_type,
            "partition": self.partition,
            "generated_at": datetime.now().isoformat(),
            "baseline": {
                "nodes": 1,
                "time_seconds": self._baseline_time or 0.0
            },
            "results": [
                {
                    "nodes": r.node_count,
                    "total_ranks": r.total_mpi_ranks,
                    "time_seconds": r.wall_time_seconds,
                    "speedup": r.speedup,
                    "efficiency": r.efficiency,
                    "status": r.status
                }
                for r in sorted_results
            ]
        }
        
        # Write JSON
        filename = f"{self.scaling_type}_scaling_report.json"
        filepath = self.data_dir / filename
        with open(filepath, 'w') as f:
            json.dump(report, f, indent=2)
        
        logger.info(f"✓ Wrote aggregated report: {filepath}")
        
        # Also write text report
        self._write_text_report(sorted_results)
        
        return filepath
    
    def _write_text_report(self, sorted_results: List[RunResult]) -> Path:
        """Write human-readable text report."""
        lines = [
            "=" * 70,
            f" {self.scaling_type.upper()} SCALING REPORT",
            "=" * 70,
            f"",
            f"Test: {self.test_name}",
            f"Partition: {self.partition}",
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"",
            "-" * 70,
        ]
        
        if self.scaling_type == "strong":
            lines.extend([
                f"{'Nodes':>8} {'Ranks':>8} {'Time(s)':>12} {'Speedup':>10} {'Efficiency':>12}",
                "-" * 70,
            ])
            
            for r in sorted_results:
                lines.append(
                    f"{r.node_count:>8} {r.total_mpi_ranks:>8} "
                    f"{r.wall_time_seconds:>12.2f} {r.speedup:>10.2f} "
                    f"{r.efficiency:>11.1f}%"
                )
        else:
            lines.extend([
                f"{'Nodes':>8} {'Ranks':>8} {'Time(s)':>12} {'Efficiency':>12}",
                "-" * 70,
            ])
            
            for r in sorted_results:
                lines.append(
                    f"{r.node_count:>8} {r.total_mpi_ranks:>8} "
                    f"{r.wall_time_seconds:>12.2f} {r.efficiency:>11.1f}%"
                )
        
        lines.extend([
            "-" * 70,
            "",
            "Notes:",
            f"  - Baseline time: {self._baseline_time or 0:.2f}s (1 node)",
        ])
        
        if self.scaling_type == "strong":
            lines.append("  - Ideal speedup = N (linear)")
            lines.append("  - Efficiency = speedup / N × 100%")
        else:
            lines.append("  - Ideal efficiency = 100% (constant time)")
            lines.append("  - Efficiency = T1 / Tn × 100%")
        
        lines.append("=" * 70)
        
        # Write file
        filename = f"{self.scaling_type}_scaling_report.txt"
        filepath = self.data_dir / filename
        with open(filepath, 'w') as f:
            f.write("\n".join(lines))
        
        logger.info(f"✓ Wrote text report: {filepath}")
        
        return filepath
    
    def get_result(self, node_count: int) -> Optional[RunResult]:
        """Get result for a specific node count."""
        return self._results.get(node_count)
    
    def get_all_results(self) -> List[RunResult]:
        """Get all results sorted by node count."""
        return sorted(self._results.values(), key=lambda r: r.node_count)


# =============================================================================
# Convenience functions
# =============================================================================

def create_result_writer(
    output_dir: Path,
    scaling_type: str,
    partition: str = "",
    test_name: str = "scaling_test"
) -> IncrementalResultWriter:
    """Create an incremental result writer."""
    return IncrementalResultWriter(output_dir, scaling_type, partition, test_name)


# =============================================================================
# Self-test
# =============================================================================

if __name__ == '__main__':
    import tempfile
    
    logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
    
    print("=" * 70)
    print(" Incremental Result Writer - Self Test")
    print("=" * 70)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        writer = IncrementalResultWriter(
            output_dir=Path(tmpdir),
            scaling_type="weak",
            partition="boost_usr_prod",
            test_name="test_scaling"
        )
        
        # Simulate runs
        for nodes in [1, 2, 4]:
            result = writer.create_pending_result(
                node_count=nodes,
                total_mpi_ranks=nodes * 4,
                mpi_ranks_per_node=4,
                cores_per_rank=8,
                gpus_per_node=4,
                procs_decomp=(nodes, 2, 2),
                job_directory=f"{tmpdir}/node{nodes}",
                slurm_job_id=f"12345{nodes}"
            )
            
            # Simulate completion
            writer.update_run_status(
                node_count=nodes,
                status="completed",
                wall_time=100.0 + nodes * 5,  # Simulated time
                exit_code=0
            )
        
        # Write summary
        writer.write_summary()
        writer.write_aggregated_report()
        
        # Verify files
        data_dir = Path(tmpdir) / "data"
        expected_files = [
            "weak_scaling_report_node1.json",
            "weak_scaling_report_node2.json",
            "weak_scaling_report_node4.json",
            "summary.json",
            "weak_scaling_report.json",
            "weak_scaling_report.txt"
        ]
        
        print("\n  Checking generated files:")
        for filename in expected_files:
            filepath = data_dir / filename
            if filepath.exists():
                print(f"    ✓ {filename}")
            else:
                print(f"    ✗ {filename} MISSING")
                raise AssertionError(f"Missing file: {filename}")
        
        # Verify content
        with open(data_dir / "summary.json") as f:
            summary = json.load(f)
        
        assert summary['completed_runs'] == 3
        assert summary['total_runs'] == 3
        print(f"\n    Summary: {summary['completed_runs']}/{summary['total_runs']} runs completed")
    
    print("\n" + "=" * 70)
    print(" All tests PASSED")
    print("=" * 70)
