"""
Local scheduler backend for running jobs as local processes.
"""

import subprocess
import logging
from pathlib import Path
from typing import List, Optional
import time
import signal

from core.abstracts import SchedulerInterface
from core.config import JobConfig, ResourceConfig
from core.types import JobStatus


logger = logging.getLogger(__name__)


class LocalScheduler(SchedulerInterface):
    """Run jobs as local processes."""
    
    def __init__(self, options: Optional[dict] = None):
        super().__init__(options)
        self.processes = {}  # job_id -> subprocess.Popen
    
    def generate_job_script(
        self,
        job_config: JobConfig,
        resource_config: ResourceConfig,
        command: List[str],
        env_setup: List[str]
    ) -> str:
        """Generate a script for local execution."""
        import platform
        
        # Check if we're on Windows
        is_windows = platform.system() == "Windows"
        
        if is_windows:
            # Windows batch script
            script = "@echo off\n"
            script += "echo ======================================\n"
            script += "echo ======================================\n"
            script += "echo This is my script\n"
            script += "type %0\n"
            script += "echo ======================================\n"
            script += "echo ======================================\n\n"
        else:
            # Unix bash script
            script = "#!/bin/bash\n\n"
            script += "echo \"======================================\"\n"
            script += "echo \"======================================\"\n"
            script += "echo \"This is my script\"\n"
            script += "cat $0\n"
            script += "echo \"======================================\"\n"
            script += "echo \"======================================\"\n\n"
        
        # Environment setup
        for cmd in env_setup:
            if is_windows:
                # Convert export to set for Windows
                if cmd.startswith("export "):
                    var_part = cmd[7:]  # Remove "export "
                    if "=" in var_part:
                        var_name, var_value = var_part.split("=", 1)
                        script += f"set {var_name}={var_value}\n"
                    else:
                        script += f"set {var_part}\n"
                else:
                    script += f"{cmd}\n"
            else:
                script += f"{cmd}\n"
        
        if is_windows:
            script += '\necho Loading modules\n\n'
        else:
            script += '\necho "Loading modules"\n\n'
        
        # Change to working directory
        if job_config.working_dir:
            if is_windows:
                # Use relative path to avoid issues with absolute paths
                script += f"cd /d \"%~dp0\"\n\n"
            else:
                script += f"cd {job_config.working_dir}\n\n"
        
        if is_windows:
            script += "set OMP_NUM_THREADS=1\n"
        else:
            script += "export OMP_NUM_THREADS=1\n"
        
        # Extract binary path as variable for clarity (same as Slurm)
        command_copy = command.copy() if command else []
        
        if command_copy and len(command_copy) > 0 and not is_windows:
            # Find the executable in the command (must contain '/' and not be a flag)
            binary_idx = None
            for idx, arg in enumerate(command_copy):
                # Skip flags (anything starting with -)
                if arg.startswith('-'):
                    continue
                # Skip launcher commands
                if arg in ['srun', 'mpirun', 'mpiexec']:
                    continue
                # Must be a path (contains /) and not a number
                if '/' in arg and not arg.replace('.', '').replace('/', '').isdigit():
                    binary_idx = idx
                    break
            
            if binary_idx is not None:
                from pathlib import Path
                full_binary_path = command_copy[binary_idx]
                binary_path_obj = Path(full_binary_path)
                binary_dir = str(binary_path_obj.parent)
                binary_name = binary_path_obj.name
                
                script += f"\n# Binary location\n"
                script += f"BINARY={binary_dir}\n"
                # Replace in command with $BINARY/name
                command_copy[binary_idx] = f"$BINARY/{binary_name}"
        
        if is_windows:
            script += "echo %date% %time%\n\n"
        else:
            script += "\ndate\n\n"
        
        # Main command with timing
        if is_windows:
            script += 'set start_time=%time%\n'
            script += "echo Execute command\n"
            script += " ".join(command) + "\n\n"
        else:
            script += 'start_time="$(date -u +%s.%N)"\n'
            script += "# Execute command\n"
            script += " ".join(command_copy) + "\n\n"
        
        if is_windows:
            script += 'set end_time=%time%\n'
            script += 'echo Total time elapsed for process\n'
        else:
            script += 'end_time="$(date -u +%s.%N)"\n'
            script += 'elapsed="$(bc <<<"$end_time-$start_time")"\n'
            script += 'echo "Total of $elapsed seconds elapsed for process"\n'
        
        if is_windows:
            script += "echo %date% %time%\n"
        else:
            script += "date\n"
        
        return script
    
    def submit_job(self, script_path: Path) -> str:
        """Execute script as a local process."""
        import platform
        import os
        
        try:
            # Create unique job ID
            job_id = f"local_{int(time.time() * 1000)}"
            
            # Get output and error paths
            output_path = script_path.parent / "out_local"
            error_path = script_path.parent / "err_local"
            
            # Determine the command to run based on the OS
            is_windows = platform.system() == "Windows"
            if is_windows:
                # On Windows, use cmd.exe for batch files
                cmd = ["cmd.exe", "/c", str(script_path)]
            else:
                # On Unix-like systems, make sure the script is executable and run it directly
                script_path.chmod(0o755)
                cmd = ["/bin/bash", str(script_path)]
            
            # Start process
            with open(output_path, 'w') as outfile, open(error_path, 'w') as errfile:
                process = subprocess.Popen(
                    cmd,
                    stdout=outfile,
                    stderr=errfile,
                    cwd=script_path.parent,
                    shell=is_windows  # Use shell on Windows
                )
            
            self.processes[job_id] = process
            logger.debug(f"Started local process {process.pid} as {job_id}")
            
            return job_id
            
        except Exception as e:
            logger.error(f"Failed to start local job: {e}")
            raise
    
    def get_job_status(self, job_id: str) -> JobStatus:
        """Check if local process is still running."""
        if job_id not in self.processes:
            return JobStatus.UNKNOWN
        
        process = self.processes[job_id]
        returncode = process.poll()
        
        if returncode is None:
            return JobStatus.RUNNING
        elif returncode == 0:
            return JobStatus.COMPLETED
        else:
            return JobStatus.FAILED
    
    def cancel_job(self, job_id: str) -> bool:
        """Terminate local process."""
        if job_id not in self.processes:
            return False
        
        try:
            process = self.processes[job_id]
            process.send_signal(signal.SIGTERM)
            time.sleep(1)
            
            if process.poll() is None:
                process.kill()
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to cancel job {job_id}: {e}")
            return False
    
    def wait_for_completion(
        self,
        job_id: str,
        timeout: Optional[int] = None
    ) -> JobStatus:
        """Wait for local process to complete."""
        if job_id not in self.processes:
            return JobStatus.UNKNOWN
        
        try:
            process = self.processes[job_id]
            returncode = process.wait(timeout=timeout)
            
            if returncode == 0:
                return JobStatus.COMPLETED
            else:
                return JobStatus.FAILED
                
        except subprocess.TimeoutExpired:
            return JobStatus.TIMEOUT
        except Exception as e:
            logger.error(f"Error waiting for job {job_id}: {e}")
            return JobStatus.FAILED
