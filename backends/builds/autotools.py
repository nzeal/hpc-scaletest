"""
GNU Autotools build system backend (./configure && make && make install).
"""

import subprocess
import logging
from pathlib import Path
from typing import Optional, Dict

from core.abstracts import BuildSystemInterface


logger = logging.getLogger(__name__)


class AutotoolsBackend(BuildSystemInterface):
    """GNU Autotools build system."""
    
    def configure(
        self,
        source_dir: Path,
        build_dir: Path,
        flags: Optional[Dict[str, str]] = None
    ) -> bool:
        """Configure with ./configure."""
        try:
            build_dir.mkdir(parents=True, exist_ok=True)
            
            configure_script = source_dir / "configure"
            if not configure_script.exists():
                logger.error(f"Configure script not found: {configure_script}")
                return False
            
            cmd = [str(configure_script), f"--prefix={build_dir / 'install'}"]
            
            # Add configure flags
            if flags:
                for key, value in flags.items():
                    if value:
                        cmd.append(f"--{key}={value}")
                    else:
                        cmd.append(f"--{key}")
            
            logger.info(f"Configuring with autotools: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                cwd=build_dir,
                capture_output=True,
                text=True,
                check=True
            )
            
            logger.debug(result.stdout)
            logger.info("Configuration completed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Configuration failed: {e.stderr}")
            return False
    
    def build(
        self,
        build_dir: Path,
        parallel_jobs: int = 1
    ) -> bool:
        """Build with make."""
        try:
            cmd = ["make", f"-j{parallel_jobs}"]
            
            logger.info(f"Building with make: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                cwd=build_dir,
                capture_output=True,
                text=True,
                check=True
            )
            
            logger.debug(result.stdout)
            logger.info("Build completed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Build failed: {e.stderr}")
            return False
    
    def install(
        self,
        build_dir: Path,
        install_dir: Path
    ) -> bool:
        """Install with make install."""
        try:
            cmd = ["make", "install"]
            
            logger.info(f"Installing: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                cwd=build_dir,
                capture_output=True,
                text=True,
                check=True
            )
            
            logger.debug(result.stdout)
            logger.info("Installation completed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Installation failed: {e.stderr}")
            return False
    
    def clean(self, build_dir: Path) -> bool:
        """Clean with make clean."""
        try:
            subprocess.run(
                ["make", "clean"],
                cwd=build_dir,
                capture_output=True,
                check=True
            )
            logger.info("Clean completed successfully")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Clean failed: {e.stderr}")
            return False
