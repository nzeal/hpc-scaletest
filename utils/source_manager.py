"""
Source Manager - HARDENED VERSION

Handles source code acquisition, primarily via Git, with robust error handling,
retry logic, and system-agnostic operations.

CRITICAL FIXES:
- Replaced direct subprocess calls with a robust command execution wrapper.
- Added retry logic for network-related Git failures.
- Improved error handling for non-zero exit codes.
- Ensured cross-platform path handling using pathlib.
"""

import subprocess
import logging
import shutil
import os
import time
from pathlib import Path
from typing import Optional, Dict, List, Tuple

logger = logging.getLogger(__name__)

# --- Configuration ---
MAX_RETRIES = 3
RETRY_DELAY_SECONDS = 5
GIT_COMMAND_TIMEOUT = 300  # 5 minutes for clone/fetch operations

# --- Custom Exceptions ---
class SourceAcquisitionError(Exception):
    """Base exception for source acquisition failures."""
    pass

class GitError(SourceAcquisitionError):
    """Raised when a Git command fails."""
    def __init__(self, message, command, returncode, stderr):
        super().__init__(message)
        self.command = command
        self.returncode = returncode
        self.stderr = stderr

class NetworkError(GitError):
    """Raised for network-related Git failures (e.g., timeout, connection reset)."""
    pass

# --- Command Execution Wrapper ---
def _execute_git_command(
    cmd: List[str], 
    cwd: Path, 
    description: str,
    timeout: int = GIT_COMMAND_TIMEOUT
) -> Tuple[str, str]:
    """
    Execute a Git command with robust error handling and retries.
    
    Args:
        cmd: Git command as list of strings (e.g., ['git', 'clone', 'url'])
        cwd: Working directory
        description: Human-readable description of the operation
        timeout: Timeout in seconds
        
    Returns:
        Tuple of (stdout, stderr)
        
    Raises:
        GitError: If the command fails after all retries.
    """
    cmd_str = ' '.join(cmd)
    
    for attempt in range(MAX_RETRIES):
        logger.info(f"Attempt {attempt + 1}/{MAX_RETRIES}: Executing Git command: {description} in {cwd}")
        
        try:
            result = subprocess.run(
                cmd,
                cwd=str(cwd),
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False, # Do not raise exception on non-zero exit code yet
                env=os.environ.copy()
            )
            
            # Check for common network errors in stderr/stdout
            network_errors = [
                'timed out', 'connection reset', 'could not resolve host', 
                'network is unreachable', 'ssh: connect to host'
            ]
            is_network_error = any(err in result.stderr.lower() or err in result.stdout.lower() for err in network_errors)
            
            if result.returncode == 0:
                logger.info(f"Git command successful: {description}")
                return result.stdout, result.stderr
            
            elif is_network_error and attempt < MAX_RETRIES - 1:
                logger.warning(f"Network error detected for {description}. Retrying in {RETRY_DELAY_SECONDS}s...")
                time.sleep(RETRY_DELAY_SECONDS)
                continue
            
            elif result.returncode != 0:
                # Non-zero exit code, not a network error, or last attempt failed
                error_message = f"Git command failed (Exit Code {result.returncode}): {description}"
                logger.error(error_message)
                logger.error(f"Command: {cmd_str}")
                logger.error(f"Stderr: {result.stderr.strip()}")
                
                # Check for authentication failure
                if 'authentication failed' in result.stderr.lower() or 'permission denied' in result.stderr.lower():
                    raise GitError(f"Authentication failed for Git repository: {description}", cmd_str, result.returncode, result.stderr)
                
                raise GitError(error_message, cmd_str, result.returncode, result.stderr)
                
        except subprocess.TimeoutExpired:
            if attempt < MAX_RETRIES - 1:
                logger.warning(f"Git command timed out for {description}. Retrying in {RETRY_DELAY_SECONDS}s...")
                time.sleep(RETRY_DELAY_SECONDS)
                continue
            else:
                raise NetworkError(f"Git command timed out after {MAX_RETRIES} attempts: {description}", cmd_str, -1, "TimeoutExpired")
        
        except FileNotFoundError:
            raise GitError(f"Git executable not found. Is Git installed and in PATH?", cmd_str, -1, "Git executable not found")
        
        except Exception as e:
            # Catch all other exceptions (e.g., permission errors)
            # We only catch unexpected exceptions here, not the GitError raised above
            if isinstance(e, GitError):
                raise e
            raise GitError(f"Unexpected error during Git operation: {e}", cmd_str, -1, str(e))

    # Should be unreachable, but for safety
    raise GitError(f"Git command failed after {MAX_RETRIES} attempts: {description}", cmd_str, -1, "Unknown failure after retries")


# --- SourceManager Class ---
class SourceManager:
    """
    Manages the acquisition and preparation of source code.
    """
    
    def __init__(self):
        logger.info("SourceManager initialized.")
    
    def acquire_source(
        self, 
        repo_url: str, 
        target_dir: Path, 
        branch: Optional[str] = None,
        tag: Optional[str] = None,
        commit: Optional[str] = None,
        shallow: bool = True,
        submodules: bool = False
    ) -> Path:
        """
        Acquires source code from a Git repository.
        
        Args:
            repo_url: The URL of the Git repository.
            target_dir: The directory where the source code should be placed.
            branch: The branch to checkout.
            tag: The tag to checkout.
            commit: The specific commit hash to checkout.
            shallow: If True, performs a shallow clone (depth 1).
            submodules: If True, initializes and updates submodules.
            
        Returns:
            The path to the acquired source code directory.
            
        Raises:
            SourceAcquisitionError: If the acquisition fails.
        """
        target_dir = Path(target_dir).resolve()
        
        if target_dir.exists():
            logger.warning(f"Target directory {target_dir} already exists. Attempting to clean.")
            try:
                shutil.rmtree(target_dir)
            except Exception as e:
                raise SourceAcquisitionError(f"Failed to clean existing target directory {target_dir}: {e}")
        
        target_dir.mkdir(parents=True, exist_ok=True)
        
        # 1. Clone the repository
        clone_cmd = ["git", "clone"]
        if shallow:
            clone_cmd.extend(["--depth", "1"])
        
        # Use a temporary directory for the clone if a specific branch/tag/commit is needed
        # and we don't want to clone directly into the final target_dir
        temp_clone_dir = target_dir
        
        if branch or tag or commit:
            # If a specific ref is requested, we clone into a temp folder
            # and then move the contents to the final target_dir
            temp_clone_dir = target_dir.parent / f"{target_dir.name}_temp_clone"
            if temp_clone_dir.exists():
                shutil.rmtree(temp_clone_dir)
            temp_clone_dir.mkdir(parents=True, exist_ok=True)
            
            # Clone without checking out a specific branch/tag/commit yet
            clone_cmd.extend([repo_url, str(temp_clone_dir)])
            
            # Execute clone
            _execute_git_command(clone_cmd, target_dir.parent, f"Cloning {repo_url} into temp directory")
            
            # 2. Checkout the specific ref
            if commit:
                ref = commit
            elif tag:
                ref = tag
            elif branch:
                ref = branch
            else:
                ref = None
            
            if ref:
                checkout_cmd = ["git", "checkout", ref]
                _execute_git_command(checkout_cmd, temp_clone_dir, f"Checking out ref: {ref}")
            
            # 3. Move contents to final target_dir
            try:
                for item in temp_clone_dir.iterdir():
                    shutil.move(str(item), str(target_dir))
                shutil.rmtree(temp_clone_dir)
            except Exception as e:
                raise SourceAcquisitionError(f"Failed to move contents from temp dir to {target_dir}: {e}")
            
        else:
            # Simple clone directly into target_dir
            clone_cmd.extend([repo_url, str(target_dir)])
            _execute_git_command(clone_cmd, target_dir.parent, f"Cloning {repo_url} into {target_dir.name}")
        
        # 4. Handle submodules
        if submodules:
            logger.info("Initializing and updating Git submodules...")
            init_cmd = ["git", "submodule", "update", "--init", "--recursive"]
            try:
                _execute_git_command(init_cmd, target_dir, "Initializing submodules")
            except GitError as e:
                logger.warning(f"Submodule update failed: {e}. Attempting to continue.")
        
        logger.info(f"Source code successfully acquired at: {target_dir}")
        return target_dir

    def get_current_commit(self, source_dir: Path) -> str:
        """
        Gets the current commit hash of the repository.
        
        Args:
            source_dir: Path to the local Git repository.
            
        Returns:
            The full commit hash as a string.
            
        Raises:
            GitError: If the command fails.
        """
        source_dir = Path(source_dir).resolve()
        
        if not (source_dir / ".git").exists():
            raise GitError(f"Directory {source_dir} is not a Git repository.", "git rev-parse HEAD", -1, "Not a Git repository")
            
        cmd = ["git", "rev-parse", "HEAD"]
        stdout, _ = _execute_git_command(cmd, source_dir, "Getting current commit hash", timeout=10)
        return stdout.strip()

    def get_remote_url(self, source_dir: Path) -> str:
        """
        Gets the remote URL of the repository.
        
        Args:
            source_dir: Path to the local Git repository.
            
        Returns:
            The remote URL as a string.
            
        Raises:
            GitError: If the command fails.
        """
        source_dir = Path(source_dir).resolve()
        
        if not (source_dir / ".git").exists():
            raise GitError(f"Directory {source_dir} is not a Git repository.", "git remote get-url origin", -1, "Not a Git repository")
            
        cmd = ["git", "remote", "get-url", "origin"]
        stdout, _ = _execute_git_command(cmd, source_dir, "Getting remote URL", timeout=10)
        return stdout.strip()

    def update_source(self, source_dir: Path) -> str:
        """
        Updates the source code to the latest commit on the current branch.
        
        Args:
            source_dir: Path to the local Git repository.
            
        Returns:
            The new commit hash.
            
        Raises:
            GitError: If the update fails.
        """
        source_dir = Path(source_dir).resolve()
        
        if not (source_dir / ".git").exists():
            raise GitError(f"Directory {source_dir} is not a Git repository.", "git pull", -1, "Not a Git repository")
            
        # 1. Fetch latest changes
        fetch_cmd = ["git", "fetch", "origin"]
        _execute_git_command(fetch_cmd, source_dir, "Fetching latest changes")
        
        # 2. Rebase/Merge (using pull for simplicity, but fetch+reset is safer)
        # We'll use fetch + reset to ensure a clean state, avoiding merge conflicts
        
        # Get current branch name
        branch_cmd = ["git", "rev-parse", "--abbrev-ref", "HEAD"]
        branch, _ = _execute_git_command(branch_cmd, source_dir, "Getting current branch", timeout=10)
        branch = branch.strip()
        
        # Reset to the remote tracking branch
        reset_cmd = ["git", "reset", "--hard", f"origin/{branch}"]
        _execute_git_command(reset_cmd, source_dir, f"Hard resetting to origin/{branch}")
        
        # 3. Update submodules if they exist
        if (source_dir / ".gitmodules").exists():
            logger.info("Updating Git submodules after pull...")
            init_cmd = ["git", "submodule", "update", "--init", "--recursive"]
            try:
                _execute_git_command(init_cmd, source_dir, "Updating submodules")
            except GitError as e:
                logger.warning(f"Submodule update failed: {e}. Attempting to continue.")
        
        logger.info(f"Source code successfully updated in: {source_dir}")
        return self.get_current_commit(source_dir)

# Re-export exceptions for external use
SourceManager.SourceAcquisitionError = SourceAcquisitionError
SourceManager.GitError = GitError
SourceManager.NetworkError = NetworkError
