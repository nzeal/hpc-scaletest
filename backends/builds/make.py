"""
Standard Makefile backend.
"""

import subprocess
import logging
from pathlib import Path
from typing import Optional, Dict

from core.abstracts import BuildSystemInterface


logger = logging.getLogger(__name__)


class MakeBackend(BuildSystemInterface):
    def __init__(self, options: Optional[Dict] = None):
        """Initialize Make backend.
        
        Args:
            options: Optional dict with 'module_commands' list for environment setup
        """
        self.options = options or {}
        self.module_commands = self.options.get('module_commands', [])
    
    def configure(self, source_dir: Path, build_dir: Path, flags: Optional[Dict[str, str]] = None) -> bool:
        logger.info("Make: No configure step")
        return True
    
    def build(self, build_dir: Path, parallel_jobs: int = 1) -> bool:
        try:
            make_cmd = ["make", f"-j{parallel_jobs}"]
            logger.info(f"Building: {' '.join(make_cmd)}")
            
            # Run with modules loaded if provided
            if self.module_commands:
                logger.info(f"  Loading modules: {', '.join(self.module_commands)}")
                full_cmd = " && ".join(self.module_commands + [' '.join(make_cmd)])
                bash_cmd = f"bash -l -c '{full_cmd}'"
                result = subprocess.run(bash_cmd, shell=True, cwd=build_dir, capture_output=True, text=True, check=True)
            else:
                result = subprocess.run(make_cmd, cwd=build_dir, capture_output=True, text=True, check=True)
            
            logger.info("Build OK")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Build failed: {e.stderr}")
            return False
    
    def install(self, build_dir: Path, install_dir: Path) -> bool:
        try:
            make_cmd = ["make", "install", f"DESTDIR={install_dir}"]
            logger.info(f"Installing: {' '.join(make_cmd)}")
            
            # Run with modules loaded if provided
            if self.module_commands:
                full_cmd = " && ".join(self.module_commands + [' '.join(make_cmd)])
                bash_cmd = f"bash -l -c '{full_cmd}'"
                result = subprocess.run(bash_cmd, shell=True, cwd=build_dir, capture_output=True, text=True, check=True)
            else:
                result = subprocess.run(make_cmd, cwd=build_dir, capture_output=True, text=True, check=True)
            
            logger.info("Install OK")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Install failed: {e.stderr}")
            return False
    
    def clean(self, build_dir: Path) -> bool:
        try:
            subprocess.run(["make", "clean"], cwd=build_dir, capture_output=True, check=True)
            logger.info("Clean OK")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Clean failed: {e.stderr}")
            return False
