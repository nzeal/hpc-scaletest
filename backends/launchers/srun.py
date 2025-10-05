# backends/launchers/srun.py
import subprocess
import os
import re
from typing import List, Dict, Any
from ...core.abstracts import AbstractLauncher
import logging

logger = logging.getLogger(__name__)

class SrunLauncher(AbstractLauncher):
    def launch(self, command: List[str], nodes: int, procs_per_node: int, env: Dict[str, str] = None) -> Dict[str, Any]:
        """Launch parallel job using srun."""
        total_ntasks = nodes * procs_per_node
        cmd = [
            'srun',
            f'--nodes={nodes}',
            f'--ntasks={total_ntasks}',
            f'--ntasks-per-node={procs_per_node}',
            '--cpu-bind=cores'  # Basic binding; can be extended for GPUs
        ] + command

        full_env = dict(os.environ)
        if env:
            full_env.update(env)

        # Wrap with time to capture runtime
        timed_cmd = ['/usr/bin/time', '-f', '%e'] + cmd

        try:
            result = subprocess.run(
                timed_cmd,
                capture_output=True,
                text=True,
                env=full_env,
                check=True
            )
            runtime = self._parse_runtime(result.stderr)
            logger.info(f"srun launch completed with runtime {runtime}s")
            return {
                'runtime': runtime,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'returncode': result.returncode
            }
        except subprocess.CalledProcessError as e:
            logger.error(f"srun launch failed: {e.stderr}")
            return {
                'runtime': 0.0,
                'stdout': e.stdout,
                'stderr': e.stderr,
                'returncode': e.returncode,
                'status': 'FAILED'
            }

    def _parse_runtime(self, output: str) -> float:
        """Parse wall-clock time from /usr/bin/time output."""
        match = re.search(r'(\d+\.?\d*)', output.strip())
        if match:
            return float(match.group(1))
        logger.warning("Could not parse runtime from output")
        return 0.0
