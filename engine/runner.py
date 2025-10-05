# engine/runner.py
import os
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Any
from ..core.test_definition import Test
from ..core.factory import BackendFactory
from ..core.types import TestConfig, JobID, JobStatus
from ..utils.file_utils import modify_file
from .scaling import ScalingEngine
import subprocess
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

class TestRunner:
    def __init__(self, test: Test):
        self.test = test
        self.config: TestConfig = test.get_config()
        self.scheduler = BackendFactory.create_scheduler(self.config.scheduler)
        self.launcher = BackendFactory.create_launcher(self.config.launcher)
        self.module_system = BackendFactory.create_module_system(self.config.module_system)
        self.build_system = None
        if self.config.build_system:
            self.build_system = BackendFactory.create_build_system(self.config.build_system)
        self.scaling_engine = ScalingEngine()
        self.template_jobscript = 'job.sh'

    def run(self, output_dir: Path = None):
        """Run the scaling test."""
        if output_dir is None:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            output_dir = Path(f"output/{self.config.name}_{self.config.scaling.scaling_type.name.lower()}_{timestamp}")
        output_dir.mkdir(parents=True, exist_ok=True)

        # Build if needed
        if self.build_system:
            logger.info("Building application...")
            source_dir = Path("./src")  # Assume src dir
            build_path = self.build_system.build(source_dir)
            install_path = self.build_system.install(build_path)
            # Update command to point to installed exe
            self.config.command[0] = str(install_path / self.config.command[0])

        # Setup modules and env globally if not per-job
        self.module_system.purge()
        self.module_system.load(self.config.modules)
        # Env is set in job script

        # Generate scaling configs
        scaling_configs = self.scaling_engine.generate_scaling_configs(
            self.config.scaling, self.config.resources.procs_per_node
        )

        results = {}
        for actual_nodes, sc in scaling_configs.items():
            node_dir = output_dir / f"nodes_{actual_nodes}"
            node_dir.mkdir(exist_ok=True)

            # Copy input file to node dir
            input_copy = node_dir / self.config.input_file.name
            self.config.input_file.copy(input_copy)

            # Generate job script
            job_script = self._generate_job_script(node_dir, sc)

            # Modify for large jobs
            if actual_nodes > 64:
                logger.info(f"Applying boost QoS for large job: {actual_nodes} nodes")
                modify_file(
                    job_script,
                    [
                        ('###SBATCH --qos=boost_qos_bprod', '#SBATCH --qos=boost_qos_bprod')
                    ]
                )

            # Submit and wait for completion
            logger.info(f"Submitting job for {actual_nodes} nodes...")
            job_id: JobID = self.scheduler.submit(job_script)
            job_result = self.scheduler.wait(job_id, timeout=3600)  # 1 hour default

            # Parse results from job_result
            runtime = self._parse_runtime_from_result(job_result)
            results[actual_nodes] = {
                'runtime': runtime,
                'status': job_result.get('status', JobStatus.FAILED),
                'total_procs': sc['total_procs']
            }

            # Save per-node results
            with open(node_dir / 'results.json', 'w') as f:
                json.dump(results[actual_nodes], f, indent=2)

        # Aggregate and analyze
        self._aggregate_results(output_dir, results, scaling_configs)
        logger.info(f"Test run completed. Results in {output_dir}")

    def _generate_job_script(self, node_dir: Path, sc: Dict) -> Path:
        """Generate the job script for a specific scaling config."""
        job_script = node_dir / self.template_jobscript
        content = self._get_job_template(node_dir, sc)
        job_script.write_text(content)
        os.chmod(job_script, 0o755)
        return job_script

    def _get_job_template(self, node_dir: Path, sc: Dict) -> str:
        """Generate content for job.sh based on backend."""
        n_nodes = sc['n_nodes']
        gpus_line = f"#SBATCH --gres=gpu:{self.config.resources.gpus_per_node}" if self.config.resources.gpus_per_node > 0 else ""
        env_vars = "\n".join([f'export {k}="{v}"' for k, v in self.config.env.items()])
        module_loads = "\n".join([f'module load {m}' for m in self.config.modules])

        # Generate launcher command
        if isinstance(self.launcher, SrunLauncher):
            launch_cmd = f"srun --nodes={n_nodes} --ntasks-per-node={self.config.resources.procs_per_node} --cpu-bind=cores"
        else:  # MpiRun
            launch_cmd = f"mpirun -np {sc['total_procs']}"

        # Benchmark command args
        input_arg = f"--input {self.config.input_file.name}"
        scaling_args = " ".join([
            f"--procs-x {sc['xproc']}",
            f"--procs-y {sc['yproc']}",
            f"--procs-z {sc['zproc']}",
            f"--domain-x {sc['lx']}",
            f"--domain-y {sc['ly']}",
            f"--domain-z {sc['lz']}",
            f"--cells-x {sc['nxc']}",
            f"--cells-y {sc['nyc']}",
            f"--cells-z {sc['nzc']}"
        ])
        full_bench_cmd = f"{self.config.command[0]} {input_arg} {scaling_args}"

        # Wrap with time for runtime
        timed_cmd = f"/usr/bin/time -f 'Runtime: %e s' {launch_cmd} {full_bench_cmd}"

        template = f"""#!/bin/bash
#SBATCH --job-name={self.config.name}_{n_nodes}n
#SBATCH --nodes={n_nodes}
#SBATCH --ntasks-per-node={self.config.resources.procs_per_node}
{gpus_line}
#SBATCH --mem={self.config.resources.memory_per_node}
#SBATCH --time={self.config.resources.time_limit}
#SBATCH --partition={self.config.resources.partition}
#SBATCH --account={self.config.resources.account}
#SBATCH --qos=normal
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

{env_vars}

{module_loads}

cd {node_dir}

# Scaling parameters
export XPROC={sc['xproc']}
export YPROC={sc['yproc']}
export ZPROC={sc['zproc']}
export LX={sc['lx']}
export LY={sc['ly']}
export LZ={sc['lz']}
export NXC={sc['nxc']}
export NYC={sc['nyc']}
export NZC={sc['nzc']}

{timed_cmd}
echo "Job completed for {n_nodes} nodes"
"""
        return template

    def _parse_runtime_from_result(self, job_result: Dict[str, Any]) -> float:
        """Parse runtime from job result stdout."""
        if 'stdout' in job_result:
            import re
            match = re.search(r'Runtime:\s*(\d+\.?\d*)\s*s', job_result['stdout'])
            if match:
                return float(match.group(1))
        elif 'parsed_runtime' in job_result:
            return job_result['parsed_runtime']
        logger.warning("Could not parse runtime; defaulting to 0")
        return 0.0

    def _aggregate_results(self, output_dir: Path, results: Dict[int, Dict], scaling_configs: Dict):
        """Aggregate results, compute efficiency, generate report and plot."""
        summary = {'results': results, 'scaling_type': self.config.scaling.scaling_type.name}
        
        if self.config.scaling.scaling_type == ScalingType.STRONG:
            initial_runtime = results.get(min(results.keys()), {}).get('runtime', 1.0)
            initial_procs = self.config.scaling.initial_procs[0] * self.config.scaling.initial_procs[1] * self.config.scaling.initial_procs[2]
            for n_nodes, res in results.items():
                runtime = res['runtime']
                total_procs = res['total_procs']
                speedup = initial_runtime / runtime if runtime > 0 else 0
                efficiency = (speedup / (total_procs / initial_procs)) * 100 if initial_procs > 0 else 0
                res['speedup'] = speedup
                res['efficiency'] = efficiency

        # Save summary
        with open(output_dir / 'summary.json', 'w') as f:
            json.dump(summary, f, indent=2)

        # Generate text report
        self._generate_efficiency_report(output_dir, results)

        # Generate plot
        self._plot_efficiency(output_dir, results)

    def _generate_efficiency_report(self, output_dir: Path, results: Dict[int, Dict]):
        """Generate text-based efficiency report."""
        report_path = output_dir / 'efficiency_report.txt'
        with open(report_path, 'w') as f:
            f.write(f"{self.config.scaling.scaling_type.name} Scaling Efficiency Report\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"{'Nodes':<8} {'Procs':<8} {'Time(s)':<8} {'Speedup':<8} {'Efficiency':<10}\n")
            f.write("-" * 50 + "\n")
            initial_procs = self.config.scaling.initial_procs[0] * self.config.scaling.initial_procs[1] * self.config.scaling.initial_procs[2]
            initial_time = results.get(1, {}).get('runtime', 1.0)
            for n_nodes in sorted(results.keys()):
                res = results[n_nodes]
                total_procs = res['total_procs']
                runtime = res['runtime']
                speedup = initial_time / runtime if runtime > 0 else 0
                efficiency = (speedup / (total_procs / initial_procs)) * 100 if total_procs > 0 else 0
                f.write(f"{n_nodes:<8} {total_procs:<8} {runtime:<8.2f} {speedup:<8.2f} {efficiency:<10.1f}%\n")

    def _plot_efficiency(self, output_dir: Path, results: Dict[int, Dict]):
        """Generate efficiency plot using matplotlib."""
        if self.config.scaling.scaling_type != ScalingType.STRONG:
            logger.info("Skipping efficiency plot for weak scaling")
            return
        nodes = sorted(results.keys())
        runtimes = [results[n]['runtime'] for n in nodes]
        plt.figure(figsize=(10, 6))
        plt.subplot(1, 2, 1)
        plt.plot(nodes, runtimes, marker='o')
        plt.xlabel('Nodes')
        plt.ylabel('Runtime (s)')
        plt.title('Runtime vs Nodes')
        plt.loglog()  # Log scale for scaling
        plt.subplot(1, 2, 2)
        efficiencies = [results[n]['efficiency'] for n in nodes]
        plt.plot(nodes, efficiencies, marker='o', color='green')
        plt.xlabel('Nodes')
        plt.ylabel('Efficiency (%)')
        plt.title('Efficiency vs Nodes')
        plt.savefig(output_dir / 'scaling_plots.png', dpi=150, bbox_inches='tight')
        plt.close()
        logger.info("Generated scaling plots")
