# backends/modules/tmod4.py
import subprocess
import os
from typing import List
from ...core.abstracts import AbstractModuleSystem
import logging

logger = logging.getLogger(__name__)

class TMod4Backend(AbstractModuleSystem):
    """Backend for TCL-based Environment Modules (version 4.x), with improved error handling."""

    def load(self, modules: List[str]) -> None:
        """Load modules using 'module load'; supports version 4.x syntax."""
        env = dict(os.environ)
        for module in modules:
            try:
                result = subprocess.run(
                    ['module', 'load', module],
                    env=env,
                    check=True,
                    capture_output=True,
                    text=True
                )
                logger.info(f"TMod4: Loaded module {module}")
                env.update(os.environ)
            except subprocess.CalledProcessError as e:
                logger.error(f"TMod4: Failed to load module {module}: {e.stderr}")
                raise

    def unload(self, modules: List[str]) -> None:
        """Unload modules using 'module unload'."""
        env = dict(os.environ)
        for module in modules:
            try:
                result = subprocess.run(
                    ['module', 'unload', module],
                    env=env,
                    check=True,
                    capture_output=True,
                    text=True
                )
                logger.info(f"TMod4: Unloaded module {module}")
                env.update(os.environ)
            except subprocess.CalledProcessError as e:
                logger.error(f"TMod4: Failed to unload module {module}: {e.stderr}")
                raise

    def purge(self) -> None:
        """Purge all modules using 'module purge'."""
        try:
            result = subprocess.run(
                ['module', 'purge'],
                check=True,
                capture_output=True,
                text=True
            )
            logger.info("TMod4: Purged all modules")
        except subprocess.CalledProcessError as e:
            logger.error(f"TMod4: Failed to purge modules: {e.stderr}")
            raise

    def list(self) -> List[str]:
        """List currently loaded modules using 'module list -t'."""
        try:
            result = subprocess.run(
                ['module', 'list', '-t'],
                capture_output=True,
                text=True,
                check=True
            )
            modules = [line.strip() for line in result.stdout.splitlines() if line.strip() and not line.startswith('Currently')]
            logger.debug(f"TMod4: Listed {len(modules)} modules")
            return modules
        except subprocess.CalledProcessError as e:
            logger.error(f"TMod4: Failed to list modules: {e.stderr}")
            return []
