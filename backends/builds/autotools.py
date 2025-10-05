# backends/builds/autotools.py
import subprocess
import os
from pathlib import Path
from typing import Dict
from ...core.abstracts import AbstractBuildSystem
import logging

logger = logging.getLogger(__name__)

class AutotoolsBackend(AbstractBuildSystem):
    def build(self, source_dir: Path, flags: Dict[str, str] = None) -> Path:
        """Build using autotools (configure/make)."""
        build_dir = source_dir / 'build'
        build_dir.mkdir(exist_ok=True)
        os.chdir(build_dir)
        configure_cmd = ['../configure']
        env = dict(os.environ)
        if flags:
            for key, value in flags.items():
                configure_cmd.extend([f'--{key}={value}'])
                logger.info(f"Setting configure flag --{key}={value}")
        try:
            subprocess.run(configure_cmd, env=env, check=True, capture_output=True, text=True)
            make_cmd = ['make', '-j$(nproc)']
            subprocess.run(make_cmd, env=env, check=True, capture_output=True, text=True)
            logger.info("Autotools build successful")
            return build_dir
        except subprocess.CalledProcessError as e:
            logger.error(f"Autotools build failed: {e.stderr}")
            raise

    def install(self, build_path: Path) -> Path:
        """Install using make install."""
        os.chdir(build_path)
        install_dir = os.environ.get('PREFIX', '/usr/local')
        cmd = ['make', 'install', f'prefix={install_dir}']
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"Autotools install successful to {install_dir}")
            return Path(install_dir)
        except subprocess.CalledProcessError as e:
            logger.error(f"Autotools install failed: {e.stderr}")
            raise
