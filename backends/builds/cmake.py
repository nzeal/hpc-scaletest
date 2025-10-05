# backends/builds/cmake.py
import subprocess
import os
from pathlib import Path
from typing import Dict
from ...core.abstracts import AbstractBuildSystem
import logging

logger = logging.getLogger(__name__)

class CMakeBackend(AbstractBuildSystem):
    def build(self, source_dir: Path, flags: Dict[str, str] = None) -> Path:
        """Build using CMake."""
        build_dir = source_dir / 'build'
        build_dir.mkdir(exist_ok=True)
        os.chdir(build_dir)
        cmake_cmd = ['cmake', str(source_dir)]
        if flags:
            for key, value in flags.items():
                cmake_cmd.extend([f'-D{key}={value}'])
                logger.info(f"Setting CMake flag {key}={value}")
        try:
            subprocess.run(cmake_cmd, check=True, capture_output=True, text=True)
            make_cmd = ['make', '-j$(nproc)']
            subprocess.run(make_cmd, check=True, capture_output=True, text=True)
            logger.info("CMake build successful")
            return build_dir
        except subprocess.CalledProcessError as e:
            logger.error(f"CMake build failed: {e.stderr}")
            raise

    def install(self, build_path: Path) -> Path:
        """Install using make install."""
        os.chdir(build_path)
        install_dir = os.environ.get('CMAKE_INSTALL_PREFIX', '/usr/local')
        cmd = ['make', 'install', f'DESTDIR={install_dir}']
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"CMake install successful to {install_dir}")
            return Path(install_dir)
        except subprocess.CalledProcessError as e:
            logger.error(f"CMake install failed: {e.stderr}")
            raise
