# backends/modules/lmod.py
import subprocess
import os
from typing import List
from ...core.abstracts import AbstractModuleSystem
import logging

logger = logging.getLogger(__name__)

class LModBackend(AbstractModuleSystem):
    """Backend for Lua-based Lmod module system."""

    def load(self, modules: List[str]) -> None:
        """Load modules using 'module load'."""
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
                logger.info(f"LMod: Loaded module {module}")
                env.update(os.environ)
            except subprocess.CalledProcessError as e:
                logger.error(f"LMod: Failed to load module {module}: {e.stderr}")
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
                logger.info(f"LMod: Unloaded module {module}")
                env.update(os.environ)
            except subprocess.CalledProcessError as e:
                logger.error(f"LMod: Failed to unload module {module}: {e.stderr}")
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
            logger.info("LMod: Purged all modules")
        except subprocess.CalledProcessError as e:
            logger.error(f"LMod: Failed to purge modules: {e.stderr}")
            raise

    def list(self) -> List[str]:
        """List currently loaded modules using 'module list'."""
        try:
            result = subprocess.run(
                ['module', 'list'],
                capture_output=True,
                text=True,
                check=True
            )
            # Lmod output format: parse lines after header
            modules = []
            in_list = False
            for line in result.stdout.splitlines():
                if 'Currently' in line:
                    in_list = True
                    continue
                if in_list and line.strip():
                    modules.append(line.strip().split()[0])  # First word is module name
            logger.debug(f"LMod: Listed {len(modules)} modules")
            return modules
        except subprocess.CalledProcessError as e:
            logger.error(f"LMod: Failed to list modules: {e.stderr}")
            return []

    def spider(self, query: str) -> List[str]:
        """Search for modules using 'module spider' (Lmod-specific)."""
        try:
            result = subprocess.run(
                ['module', 'spider', query],
                capture_output=True,
                text=True,
                check=True
            )
            # Parse spider output for matching modules
            modules = []
            for line in result.stdout.splitlines():
                if query.lower() in line.lower() and 'No' not in line:
                    modules.append(line.strip())
            logger.info(f"LMod: Spider search for '{query}' found {len(modules)} results")
            return modules
        except subprocess.CalledProcessError as e:
            logger.error(f"LMod: Spider search failed: {e.stderr}")
            return []
