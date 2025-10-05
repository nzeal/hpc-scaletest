# backends/launchers/mpirun.py
import subprocess
import os
import re
from typing import List, Dict, Any
from ...core.abstracts import AbstractLauncher
import logging

logger = logging.getLogger(__name__)

class MpiRunLauncher(AbstractLauncher):
    def launch(self, command: List[str], nodes: int, procs_per_node: int, env: Dict[str, str] = None) -> Dict[str, Any]:
        """Launch parallel job using mpirun."""
        total_procs = nodes * procs_per_node
        cmd = ['mpirun', f'-np', str(total_procs), '-hostfile', '/dev/null'] + command  # Basic; hostfile can be generated if needed

        full_env = dict(os.environ)
        if env:
            full_env.update(env)

        # Wrap with /usr/bin/time -f "real %e" for timing
        timed_cmd = ['/usr/bin/time', '-f', 'real %e'] + cmd

        try:
            result = subprocess.run(
                timed_cmd,
                capture_output=True,
                text=True,
                env=full_env,
                check=True
            )
            runtime = self._parse_runtime(result.stderr)
            logger.info(f"mpirun launch completed with runtime {runtime}s")
            return {
                'runtime': runtime,
                'stdout': result.stdout,
                'stderr': result.stderr,
                'returncode': result.returncode
            }
        except subprocess.CalledProcessError as e:
            logger.error(f"mpirun launch failed: {e.stderr}")
            return {
                'runtime': 0.0,
                'stdout': e.stdout,
                'stderr': e.stderr,
                'returncode': e.returncode,
                'status': 'FAILED'
            }

    def _parse_runtime(self, output: str) -> float:
        """Parse wall-clock time from /usr/bin/time output."""
        match = re.search(r'real\s+(\d+\.\d+)', output)
        if match:
            return float(match.group(1))
        logger.warning("Could not parse runtime from output")
        return 0.0
