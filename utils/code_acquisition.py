"""
Code acquisition utilities for fetching HPC application source code.
Supports local paths and Git repositories.
"""

import logging
import subprocess
import shutil
from pathlib import Path
from typing import Optional, Tuple
from urllib.parse import urlparse

logger = logging.getLogger(__name__)


class CodeAcquisition:
    """Handles fetching code from local paths or Git repositories."""
    
    def __init__(self, workspace_dir: Path = Path("workspace")):
        """
        Initialize code acquisition handler.
        
        Args:
            workspace_dir: Directory where Git repos will be cloned
        """
        self.workspace_dir = workspace_dir
        self.workspace_dir.mkdir(parents=True, exist_ok=True)
    
    def acquire(self, source: str) -> Tuple[Path, bool]:
        """
        Acquire source code from local path or Git repository.
        
        Args:
            source: Local path or Git repository URL
            
        Returns:
            Tuple of (source_directory, is_cloned)
            - source_directory: Path to the source code (absolute path)
            - is_cloned: True if code was cloned from Git, False if local
            
        Raises:
            ValueError: If source is invalid or inaccessible
        """
        # Check if it's a local path
        local_path = Path(source)
        if local_path.exists():
            absolute_path = local_path.resolve()
            logger.info(f"Using local source code: {absolute_path}")
            return absolute_path, False
        
        # Check if it's a Git URL
        if self._is_git_url(source):
            cloned_path = self._clone_repository(source)
            return cloned_path.resolve(), True
        
        raise ValueError(f"Invalid source: {source}. Must be a valid local path or Git URL.")
    
    def _is_git_url(self, source: str) -> bool:
        """
        Check if source is a Git URL.
        
        Args:
            source: String to check
            
        Returns:
            True if it's a Git URL
        """
        # Check common Git URL patterns
        git_patterns = [
            source.startswith("git@"),
            source.startswith("http://"),
            source.startswith("https://"),
            source.startswith("ssh://"),
            source.endswith(".git")
        ]
        return any(git_patterns)
    
    def _clone_repository(self, git_url: str) -> Path:
        """
        Clone a Git repository.
        
        Args:
            git_url: Git repository URL
            
        Returns:
            Path to cloned repository
            
        Raises:
            RuntimeError: If git clone fails
        """
        # Extract repo name from URL
        repo_name = self._extract_repo_name(git_url)
        target_dir = self.workspace_dir / repo_name
        
        # If directory already exists, remove it or use it
        if target_dir.exists():
            logger.warning(f"Repository directory already exists: {target_dir}")
            logger.info("Using existing repository. Delete it manually if you want a fresh clone.")
            return target_dir
        
        logger.info(f"Cloning repository from {git_url} to {target_dir}")
        try:
            result = subprocess.run(
                ["git", "clone", git_url, str(target_dir)],
                capture_output=True,
                text=True,
                check=True
            )
            logger.info(f"Successfully cloned repository to {target_dir}")
            logger.debug(f"Git output: {result.stdout}")
            return target_dir
        except subprocess.CalledProcessError as e:
            error_msg = f"Failed to clone repository: {e.stderr}"
            logger.error(error_msg)
            raise RuntimeError(error_msg) from e
        except FileNotFoundError:
            error_msg = "Git is not installed or not in PATH"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
    
    def _extract_repo_name(self, git_url: str) -> str:
        """
        Extract repository name from Git URL.
        
        Args:
            git_url: Git repository URL
            
        Returns:
            Repository name
        """
        # Handle different Git URL formats
        if git_url.startswith("git@"):
            # git@github.com:user/repo.git
            name = git_url.split(":")[-1]
        else:
            # https://github.com/user/repo.git
            parsed = urlparse(git_url)
            name = parsed.path.strip("/")
        
        # Remove .git extension if present
        if name.endswith(".git"):
            name = name[:-4]
        
        # Get the last component (repo name)
        name = name.split("/")[-1]
        
        return name
    
    def cleanup(self, source_dir: Path, is_cloned: bool, force: bool = False):
        """
        Clean up acquired source code.
        
        Args:
            source_dir: Path to source directory
            is_cloned: Whether the code was cloned
            force: Force cleanup even for local paths (use with caution!)
        """
        if is_cloned or force:
            if source_dir.exists():
                logger.info(f"Cleaning up source directory: {source_dir}")
                shutil.rmtree(source_dir)
        else:
            logger.info("Keeping local source directory (not cloned)")
