"""Environment module system backends."""

import subprocess
import logging
from typing import List

from core.abstracts import AbstractModuleSystem

logger = logging.getLogger(__name__)


class NoModBackend(AbstractModuleSystem):
    """No-op module backend for systems without modules."""
    
    def load(self, modules: List[str]) -> bool:
        """No-op load."""
        logger.info(f"NoModBackend: Would load modules {modules}")
        return True
    
    def unload(self, modules: List[str]) -> bool:
        """No-op unload."""
        logger.info(f"NoModBackend: Would unload modules {modules}")
        return True
    
    def purge(self) -> bool:
        """No-op purge."""
        logger.info("NoModBackend: Would purge modules")
        return True
    
    def list(self) -> List[str]:
        """No-op list."""
        return []
    
    def get_load_commands(self, modules: List[str]) -> List[str]:
        """Return empty commands (no modules to load)."""
        return []


class TModBackend(AbstractModuleSystem):
    """TCL-based module system (Environment Modules 3.x)."""
    
    def __init__(self, **options):
        self.module_cmd = options.get('module_cmd', 'modulecmd')
    
    def load(self, modules: List[str]) -> bool:
        """Load modules using module command."""
        for module in modules:
            try:
                subprocess.run(
                    [self.module_cmd, 'bash', 'load', module],
                    check=True,
                    capture_output=True
                )
                logger.info(f"Loaded module: {module}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to load module {module}: {e}")
                return False
        return True
    
    def unload(self, modules: List[str]) -> bool:
        """Unload modules."""
        for module in modules:
            try:
                subprocess.run(
                    [self.module_cmd, 'bash', 'unload', module],
                    check=True,
                    capture_output=True
                )
                logger.info(f"Unloaded module: {module}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to unload module {module}: {e}")
                return False
        return True
    
    def purge(self) -> bool:
        """Purge all modules."""
        try:
            subprocess.run(
                [self.module_cmd, 'bash', 'purge'],
                check=True,
                capture_output=True
            )
            logger.info("Purged all modules")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to purge modules: {e}")
            return False
    
    def list(self) -> List[str]:
        """List loaded modules."""
        try:
            result = subprocess.run(
                [self.module_cmd, 'bash', 'list'],
                capture_output=True,
                text=True
            )
            # Parse module list from stderr (typical for module command)
            return [line.strip() for line in result.stderr.split('\n') if line.strip()]
        except Exception as e:
            logger.error(f"Failed to list modules: {e}")
            return []
    
    def get_load_commands(self, modules: List[str]) -> List[str]:
        """Get shell commands to load modules."""
        return [f"module load {module}" for module in modules]


class TMod4Backend(TModBackend):
    """TCL-based module system (Environment Modules 4.x)."""
    
    def __init__(self, **options):
        super().__init__(**options)
        self.module_cmd = options.get('module_cmd', 'module')


class LModBackend(AbstractModuleSystem):
    """Lua-based Lmod module system."""
    
    def __init__(self, **options):
        self.module_cmd = options.get('module_cmd', 'module')
    
    def load(self, modules: List[str]) -> bool:
        """Load modules using Lmod."""
        for module in modules:
            try:
                subprocess.run(
                    ['bash', '-c', f'{self.module_cmd} load {module}'],
                    check=True,
                    capture_output=True
                )
                logger.info(f"Loaded module: {module}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to load module {module}: {e}")
                return False
        return True
    
    def unload(self, modules: List[str]) -> bool:
        """Unload modules."""
        for module in modules:
            try:
                subprocess.run(
                    ['bash', '-c', f'{self.module_cmd} unload {module}'],
                    check=True,
                    capture_output=True
                )
                logger.info(f"Unloaded module: {module}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to unload module {module}: {e}")
                return False
        return True
    
    def purge(self) -> bool:
        """Purge all modules."""
        try:
            subprocess.run(
                ['bash', '-c', f'{self.module_cmd} purge'],
                check=True,
                capture_output=True
            )
            logger.info("Purged all modules")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to purge modules: {e}")
            return False
    
    def list(self) -> List[str]:
        """List loaded modules."""
        try:
            result = subprocess.run(
                ['bash', '-c', f'{self.module_cmd} list'],
                capture_output=True,
                text=True
            )
            # Parse module list
            modules = []
            for line in result.stderr.split('\n'):
                line = line.strip()
                if line and not line.startswith('Currently') and not line.startswith('No'):
                    modules.append(line)
            return modules
        except Exception as e:
            logger.error(f"Failed to list modules: {e}")
            return []
    
    def get_load_commands(self, modules: List[str]) -> List[str]:
        """Get shell commands to load modules."""
        return [f"module load {module}" for module in modules]
    
    def spider(self, module: str) -> str:
        """Search for module in Lmod (Lmod-specific feature)."""
        try:
            result = subprocess.run(
                ['bash', '-c', f'{self.module_cmd} spider {module}'],
                capture_output=True,
                text=True
            )
            return result.stdout
        except Exception as e:
            logger.error(f"Failed to spider module {module}: {e}")
            return ""
