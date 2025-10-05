# backends/builds/spack.py
import subprocess
import os
from pathlib import Path
from typing import Dict
from ...core.abstracts import AbstractBuildSystem
import logging

logger = logging.getLogger(__name__)

class SpackBackend(AbstractBuildSystem):
    def build(self, source_dir: Path, flags: Dict[str, str] = None) -> Path:
        """Build using Spack with spec."""
        spec = flags.get('spec', 'zlib')  # Default spec if not provided
        cmd = ['spack', 'install', spec, '--jobs', '$(nproc)']
        env = dict(os.environ)
        if flags:
            for key, value in flags.items():
                if key != 'spec':
                    env[f'SPACK_{key.upper()}'] = value
                    logger.info(f"Setting Spack env SPACK_{key.upper()}={value}")
        try:
            subprocess.run(cmd, env=env, check=True, capture_output=True, text=True)
            logger.info(f"Spack install successful for {spec}")
            # Get install path using spack find --path
            path_cmd = ['spack', 'find', '--path', spec]
            path_result = subprocess.run(path_cmd, env=env, capture_output=True, text=True, check=True)
            install_path = Path(path_result.stdout.strip())
            if not install_path.exists():
                raise ValueError(f"Install path not found for {spec}: {install_path}")
            return install_path
        except subprocess.CalledProcessError as e:
            logger.error(f"Spack build failed: {e.stderr}")
            raise

    def install(self, build_path: Path) -> Path:
        """Spack combines build and install; return the build path as install."""
        logger.info(f"Spack install path: {build_path}")
        return build_path
