# backends/builds/make.py
import subprocess
import os
from pathlib import Path
from typing import Dict
from ...core.abstracts import AbstractBuildSystem
import logging

logger = logging.getLogger(__name__)

class MakeBackend(AbstractBuildSystem):
    def build(self, source_dir: Path, flags: Dict[str, str] = None) -> Path:
        """Build using make in the source directory."""
        os.chdir(source_dir)
        cmd = ['make', '-j$(nproc)']  # Use all cores
        env = dict(os.environ)
        if flags:
            for key, value in flags.items():
                env[key] = value
                logger.info(f"Setting env {key}={value} for make")
        try:
            result = subprocess.run(cmd, env=env, check=True, capture_output=True, text=True)
            logger.info("Make build successful")
            return source_dir  # Build artifacts in source_dir
        except subprocess.CalledProcessError as e:
            logger.error(f"Make build failed: {e.stderr}")
            raise

    def install(self, build_path: Path) -> Path:
        """Install using make install."""
        os.chdir(build_path)
        install_dir = os.environ.get('PREFIX', '/usr/local')
        cmd = ['make', 'install', f'PREFIX={install_dir}']
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"Make install successful to {install_dir}")
            return Path(install_dir)
        except subprocess.CalledProcessError as e:
            logger.error(f"Make install failed: {e.stderr}")
            raise
