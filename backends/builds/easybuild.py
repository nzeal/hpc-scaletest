# backends/builds/easybuild.py
import subprocess
import os
import re
from pathlib import Path
from typing import Dict
from ...core.abstracts import AbstractBuildSystem
import logging

logger = logging.getLogger(__name__)

class EasyBuildBackend(AbstractBuildSystem):
    def build(self, source_dir: Path, flags: Dict[str, str] = None) -> Path:
        """Build using EasyBuild with easyconfig file."""
        easyconfig_path = flags.get('easyconfig', str(source_dir / 'easyconfig.eb'))
        if not Path(easyconfig_path).exists():
            raise ValueError(f"Easyconfig file not found: {easyconfig_path}")
        cmd = ['eb', easyconfig_path, '--robot=1', '--job-cores=$(nproc)']
        env = dict(os.environ)
        # Pass additional flags as env if needed
        if flags:
            for key, value in flags.items():
                if key != 'easyconfig':
                    env[f'EB_{key.upper()}'] = value
                    logger.info(f"Setting EasyBuild env EB_{key.upper()}={value}")
        try:
            result = subprocess.run(cmd, env=env, check=True, capture_output=True, text=True)
            logger.info("EasyBuild build successful")
            # Parse install path from output: look for "installed into <path>"
            install_match = re.search(r'installed into\s+(\S+)', result.stdout)
            if install_match:
                build_path = Path(install_match.group(1))
            else:
                # Fallback to config
                fallback_cmd = ['eb', '--show-config', 'buildpath']
                fallback_result = subprocess.run(fallback_cmd, env=env, capture_output=True, text=True, check=True)
                build_path = Path(fallback_result.stdout.strip())
            return build_path
        except subprocess.CalledProcessError as e:
            logger.error(f"EasyBuild build failed: {e.stderr}")
            raise

    def install(self, build_path: Path) -> Path:
        """EasyBuild combines build and install; return the build path as install."""
        logger.info(f"EasyBuild install path: {build_path}")
        return build_path
