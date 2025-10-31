"""
EasyBuild backend.
"""

import subprocess
import logging
from pathlib import Path
from typing import Optional, Dict
import shutil

from core.abstracts import BuildSystemInterface


logger = logging.getLogger(__name__)


class EasyBuildBackend(BuildSystemInterface):
    def configure(self, source_dir: Path, build_dir: Path, flags: Optional[Dict[str, str]] = None) -> bool:
        logger.info("EasyBuild: Auto-config")
        return True
    
    def build(self, build_dir: Path, parallel_jobs: int = 1) -> bool:
        try:
            easyconfig = self.options.get('easyconfig', 'easyconfig.eb')
            cmd = ["eb", easyconfig, f"--parallel={parallel_jobs}", "--robot"]
            if build_dir:
                cmd += [f"--buildpath={build_dir}"]
            logger.info(f"Building: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info("Build OK")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Build failed: {e.stderr}")
            return False
    
    def install(self, build_dir: Path, install_dir: Path) -> bool:
        logger.info("EasyBuild: Auto-install")
        return True
    
    def clean(self, build_dir: Path) -> bool:
        try:
            if build_dir.exists():
                shutil.rmtree(build_dir)
            logger.info("Clean OK")
            return True
        except Exception as e:
            logger.error(f"Clean failed: {e}")
            return False
