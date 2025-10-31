"""
Spack backend.
"""

import subprocess
import logging
from pathlib import Path
from typing import Optional, Dict

from core.abstracts import BuildSystemInterface


logger = logging.getLogger(__name__)


class SpackBackend(BuildSystemInterface):
    def configure(self, source_dir: Path, build_dir: Path, flags: Optional[Dict[str, str]] = None) -> bool:
        logger.info("Spack: Auto-config")
        return True
    
    def build(self, build_dir: Path, parallel_jobs: int = 1) -> bool:
        try:
            package = self.options.get('package', 'package@version')
            cmd = ["spack", "install", f"-j{parallel_jobs}", package]
            logger.info(f"Building: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info("Build OK")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Build failed: {e.stderr}")
            return False
    
    def install(self, build_dir: Path, install_dir: Path) -> bool:
        logger.info("Spack: Auto-install")
        return True
    
    def clean(self, build_dir: Path) -> bool:
        try:
            package = self.options.get('package')
            if package:
                subprocess.run(["spack", "uninstall", "-y", package], capture_output=True, check=True)
            logger.info("Clean OK")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Clean failed: {e.stderr}")
            return False
