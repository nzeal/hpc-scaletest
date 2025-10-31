"""
CMake backend.
"""

import subprocess
import logging
from pathlib import Path
from typing import Optional, Dict

from core.abstracts import BuildSystemInterface


logger = logging.getLogger(__name__)


class CMakeBackend(BuildSystemInterface):
    def __init__(self, options: Optional[Dict] = None):
        """Initialize CMake backend.
        
        Args:
            options: Optional dict with 'module_commands' list for environment setup
        """
        self.options = options or {}
        self.module_commands = self.options.get('module_commands', [])
    def configure(self, source_dir: Path, build_dir: Path, flags: Optional[Dict[str, str]] = None) -> bool:
        try:
            # Ensure paths are absolute
            source_dir = source_dir.resolve()
            build_dir = build_dir.resolve()
            
            build_dir.mkdir(parents=True, exist_ok=True)
            cmake_cmd = ["cmake", str(source_dir)]
            if flags:
                for key, value in flags.items():
                    cmake_cmd += [f"-D{key}={value}"]
            
            logger.info(f"CMake Configure:")
            logger.info(f"  Source: {source_dir}")
            logger.info(f"  Build:  {build_dir}")
            logger.info(f"  Command: {' '.join(cmake_cmd)}")
            
            # Run with modules loaded if provided
            if self.module_commands:
                logger.info(f"  Loading modules: {', '.join(self.module_commands)}")
                full_cmd = " && ".join(self.module_commands + [' '.join(cmake_cmd)])
                bash_cmd = f"bash -l -c '{full_cmd}'"
                result = subprocess.run(bash_cmd, shell=True, cwd=build_dir, capture_output=True, text=True, check=True)
            else:
                result = subprocess.run(cmake_cmd, cwd=build_dir, capture_output=True, text=True, check=True)
            
            logger.info("Configure OK")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Configure failed:")
            logger.error(f"  Source dir: {source_dir}")
            logger.error(f"  Build dir: {build_dir}")
            logger.error(f"  Error: {e.stderr}")
            return False
    
    def build(self, build_dir: Path, parallel_jobs: int = 1) -> bool:
        try:
            cmake_cmd = ["cmake", "--build", ".", "--parallel", str(parallel_jobs)]
            logger.info(f"Building: {' '.join(cmake_cmd)}")
            
            # Run with modules loaded if provided
            if self.module_commands:
                full_cmd = " && ".join(self.module_commands + [' '.join(cmake_cmd)])
                bash_cmd = f"bash -l -c '{full_cmd}'"
                result = subprocess.run(bash_cmd, shell=True, cwd=build_dir, capture_output=True, text=True, check=True)
            else:
                result = subprocess.run(cmake_cmd, cwd=build_dir, capture_output=True, text=True, check=True)
            
            logger.info("Build OK")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Build failed: {e.stderr}")
            return False
    
    def install(self, build_dir: Path, install_dir: Path) -> bool:
        try:
            cmake_cmd = ["cmake", "--install", ".", "--prefix", str(install_dir)]
            logger.info(f"Installing: {' '.join(cmake_cmd)}")
            
            # Run with modules loaded if provided
            if self.module_commands:
                full_cmd = " && ".join(self.module_commands + [' '.join(cmake_cmd)])
                bash_cmd = f"bash -l -c '{full_cmd}'"
                result = subprocess.run(bash_cmd, shell=True, cwd=build_dir, capture_output=True, text=True, check=True)
            else:
                result = subprocess.run(cmake_cmd, cwd=build_dir, capture_output=True, text=True, check=True)
            
            logger.info("Install OK")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Install failed: {e.stderr}")
            return False
    
    def clean(self, build_dir: Path) -> bool:
        try:
            subprocess.run(["cmake", "--build", ".", "--target", "clean"], cwd=build_dir, capture_output=True, check=True)
            logger.info("Clean OK")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Clean failed: {e.stderr}")
            return False
