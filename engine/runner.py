"""
Test runner with auto-build.
"""

import logging
import json
import time
import os
import subprocess
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional

from core.test_definition import Test
from core.config import JobConfig
from core.factory import BackendFactory
from core.types import JobStatus, SchedulerBackend, LauncherBackend, ModuleBackend, BuildBackend
from core.abstracts import DefaultResultParser
from engine.scaling import ScalingEngine
from utils.file_utils import create_directory, write_file
from utils.input_analyzer import InputFileAnalyzer, validate_and_adapt_input
from utils.generic_input_parser import GenericInputParser  # NEW: Generic parser
from utils.validators import validate_job_config, ValidationError


logger = logging.getLogger(__name__)


class TestRunner:
    def __init__(self, test: Test):
        self.test = test
        test.validate()
        
        # Output dir (always use absolute path)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        scaling_type = test.scaling_config.scaling_type.value
        self.run_dir = (test.output_dir / f"{test.name}_{scaling_type}_{timestamp}").resolve()
        create_directory(self.run_dir)
        logger.info(f"Output dir: {self.run_dir}")
        
        # Track whether jobs actually completed
        self.jobs_completed = False
        
        # Load modules for build (improved for login node)
        self._load_modules_for_build()
        
        # Backends
        try:
            self.scheduler = BackendFactory.create_scheduler(
                test.backend_config.scheduler, test.backend_config.scheduler_options
            )
            logger.info(f"Created scheduler: {test.backend_config.scheduler.value}")
        except Exception as e:
            logger.error(f"Failed to create scheduler: {e}")
            raise
            
        # Determine launcher backend
        # For local scheduler, use MPIRUN by default unless SIMPLE is explicitly requested
        from core.types import LauncherBackend
        if test.backend_config.scheduler == SchedulerBackend.LOCAL:
            if test.backend_config.launcher == LauncherBackend.SIMPLE:
                launcher_backend_enum = LauncherBackend.SIMPLE
            else:
                launcher_backend_enum = LauncherBackend.MPIRUN  # Default to mpirun for local scheduler
        else:
            # Already an enum on Test.backend_config; use as-is
            launcher_backend_enum = test.backend_config.launcher
        
        logger.info(f"Selected launcher backend: {launcher_backend_enum}")
        
        try:
            self.launcher = BackendFactory.create_launcher(
                launcher_backend_enum, test.backend_config.launcher_options
            )
            logger.info(f"Created launcher: {launcher_backend_enum.value}")
        except Exception as e:
            logger.error(f"Failed to create launcher: {e}")
            raise
            
        try:
            self.module_system = BackendFactory.create_module_system(
                test.backend_config.module_system, test.backend_config.module_options
            )
            logger.info(f"Created module system: {test.backend_config.module_system.value}")
        except Exception as e:
            logger.error(f"Failed to create module system: {e}")
            raise

        # Build if source_dir set
        self.install_dir = self._build_if_needed()
        
        # Update command with built exe if built
        if self.install_dir:
            exe_name = test.build_config.executable_name
            # Convert to absolute path so job scripts work after cd
            self.test.command[0] = str((self.install_dir / exe_name).resolve())
            logger.info(f"Using binary: {self.test.command[0]}")
        
        # Create node sequence from scaling config
        node_sequence = test.scaling_config.get_node_sequence()
        
        # Extract initial values from scaling config
        initial_domain = test.scaling_config.initial_domain
        initial_cells = test.scaling_config.initial_cells
        initial_procs = test.scaling_config.initial_procs
        
        # For standalone scaling engine, we need to provide the required parameters
        # We'll create a temporary input file for the scaling engine
        self.scaling_engine = ScalingEngine(
            input_file=str(test.input_file) if test.input_file else "input.yaml",
            output_dir=str(self.run_dir / "scaling_configs"),
            nodes=node_sequence,
            procs_per_node=test.resource_config.procs_per_node,
            scale_factor=int(test.scaling_config.scaling_factor or 2.0),
            dims=test.scaling_config.scaling_dimensions,
            initial_domain=initial_domain,
            initial_cells=initial_cells,
            initial_procs=initial_procs,
            scaling_type=test.scaling_config.scaling_type.value
        )
    
    def _load_modules_for_build(self):
        """Load modules reliably on login node (bash -l for profile)."""
        modules = self.test.environment_config.modules
        if not modules:
            logger.info("No modules to load for build")
            return
        
        module_backend = self.test.backend_config.module_system
        module_system = BackendFactory.create_module_system(module_backend)
        load_cmds = module_system.generate_load_commands(modules)
        
        if load_cmds:
            full_cmd = " && ".join(load_cmds) + " && source /etc/profile"  # Source profile for env
            bash_cmd = f"bash -l -c '{full_cmd} && which mpicc'"  # Login shell, check MPI
            logger.info(f"Loading modules for build: {bash_cmd}")
            try:
                result = subprocess.run(bash_cmd, shell=True, capture_output=True, text=True, check=True)
                mpicc_path = result.stdout.strip()
                if mpicc_path and os.path.exists(mpicc_path):
                    logger.info(f"MPI found: {mpicc_path}")
                    # Set env vars for CMake
                    os.environ['CC'] = 'mpicc'
                    os.environ['CXX'] = 'mpicxx'
                    os.environ['FC'] = 'mpifort'  # If Fortran
                    os.environ['PATH'] = f"{os.path.dirname(mpicc_path)}:{os.environ.get('PATH', '')}"
                else:
                    logger.warning("MPI not found after load—build may fail (use pre-built exe)")
            except subprocess.CalledProcessError as e:
                logger.warning(f"Module load failed (non-fatal): {e.stderr}. Build may skip MPI.")
    
    def _detect_binary_in_build_dir(self, build_dir: Path) -> Optional[Path]:
        """
        Auto-detect the compiled binary in the build directory.
        Looks for executable files and returns the most likely candidate.
        """
        import stat
        
        # Common locations to check in build directory
        search_paths = [
            build_dir,                      # Root of build dir
            build_dir / "bin",              # Common bin subdirectory
            build_dir / "src",              # Source directory
            build_dir / "build",            # Sometimes nested
        ]
        
        executables = []
        
        # Also do recursive search up to 2 levels deep in main build dir
        try:
            for item in build_dir.rglob("*"):
                # Limit search depth
                if len(item.relative_to(build_dir).parts) > 2:
                    continue
                
                if item.is_file():
                    # Check if it's executable
                    try:
                        st = item.stat()
                        is_executable = bool(st.st_mode & stat.S_IXUSR)
                        
                        # Skip library files, object files, and scripts
                        skip_extensions = ['.o', '.a', '.so', '.sh', '.py', '.cmake', '.dylib']
                        if any(item.name.endswith(ext) for ext in skip_extensions):
                            continue
                        
                        # Skip common CMake/build files
                        skip_names = ['CMakeFiles', 'cmake', 'CMakeLists.txt', 'CTest', 'CPack']
                        if any(skip in item.name for skip in skip_names):
                            continue
                        
                        if is_executable and item not in executables:
                            executables.append(item)
                            logger.debug(f"Found executable: {item}")
                    except Exception as e:
                        logger.debug(f"Error checking {item}: {e}")
                        continue
        except Exception as e:
            logger.debug(f"Error during recursive search: {e}")
            
            # Fallback to simple search
            for search_path in search_paths:
                if not search_path.exists():
                    continue
                
                # Find all executable files (not .o, .a, .so, or CMake files)
                for item in search_path.iterdir():
                    if item.is_file():
                        # Check if it's executable
                        try:
                            st = item.stat()
                            is_executable = bool(st.st_mode & stat.S_IXUSR)
                            
                            # Skip library files, object files, and scripts
                            skip_extensions = ['.o', '.a', '.so', '.sh', '.py', '.cmake']
                            if any(item.name.endswith(ext) for ext in skip_extensions):
                                continue
                            
                            # Skip common CMake/build files
                            skip_names = ['CMakeFiles', 'cmake', 'CMakeLists.txt']
                            if any(skip in item.name for skip in skip_names):
                                continue
                            
                            if is_executable:
                                executables.append(item)
                                logger.debug(f"Found executable: {item}")
                        except Exception as e:
                            logger.debug(f"Error checking {item}: {e}")
                            continue
        
        if not executables:
            logger.warning("No executables found in build directory")
            return None
        
        # If multiple executables found, prefer the one that looks like the project name
        if len(executables) > 1:
            # Try to match against source directory name or test name
            source_name = self.test.build_config.source_dir.name if self.test.build_config.source_dir else ""
            test_name = self.test.name
            
            # Score each executable by how well it matches
            scored = []
            for exe in executables:
                score = 0
                exe_lower = exe.name.lower()
                
                # Prefer if it matches test name or source name
                if test_name.lower() in exe_lower:
                    score += 2
                if source_name.lower() in exe_lower:
                    score += 2
                
                # Prefer shorter names (main executable usually has simple name)
                score -= len(exe.name) * 0.01
                
                scored.append((score, exe))
            
            # Sort by score (highest first)
            scored.sort(key=lambda x: x[0], reverse=True)
            
            logger.info(f"Found {len(executables)} executables, selecting: {scored[0][1].name}")
            for score, exe in scored[:3]:  # Log top 3
                logger.debug(f"  {exe.name} (score: {score:.2f})")
            
            return scored[0][1]
        
        logger.info(f"Found executable: {executables[0].name}")
        return executables[0]
    
    def _build_if_needed(self) -> Optional[Path]:
        """Build if source_dir set, or use existing executable if available."""
        source_dir = self.test.build_config.source_dir
        if not source_dir or not source_dir.exists():
            return None
        
        # Check if build directory already exists with executable
        build_dir = self.test.build_config.build_dir or source_dir / "build"
        build_dir = build_dir.resolve()
        
        if build_dir.exists():
            logger.info(f"Build directory exists: {build_dir}")
            # Try to find existing executable
            detected_binary = self._detect_binary_in_build_dir(build_dir)
            if detected_binary and detected_binary.exists():
                logger.info(f"✓ Found existing executable: {detected_binary.name}")
                logger.info(f"  Skipping build stage - using existing binary")
                self.test.build_config.executable_name = detected_binary.name
                return build_dir
            else:
                logger.info("  No executable found in build directory, will rebuild")
        
        build_system = self.test.backend_config.build_system or BuildBackend.CMAKE
        
        # Generate module load commands to pass to build system
        build_options = dict(self.test.backend_config.build_options or {})
        modules = self.test.environment_config.modules
        if modules:
            module_backend = self.test.backend_config.module_system
            module_system = BackendFactory.create_module_system(module_backend)
            module_commands = module_system.generate_load_commands(modules)
            build_options['module_commands'] = module_commands
            logger.info(f"Build will use modules: {', '.join(modules)}")
        
        builder = BackendFactory.create_build_system(build_system, build_options)
        
        # build_dir already resolved above, just get install_dir
        install_dir = self.test.build_config.install_dir or self.run_dir / "install"
        install_dir = install_dir.resolve()  # Convert to absolute path
        
        # Flags with MPI env (now set)
        effective_flags = {
            **self.test.build_config.build_flags,
            "CMAKE_BUILD_TYPE": "Release"
        }
        
        # Only add HDF5_ROOT if environment variable is set
        if 'HDF5_ROOT' in os.environ:
            effective_flags["HDF5_ROOT"] = os.environ['HDF5_ROOT']
        elif 'HDF5_HOME' in os.environ:
            effective_flags["HDF5_ROOT"] = os.environ['HDF5_HOME']
        
        logger.info(f"Building from {source_dir} using {build_system.value}")
        if not builder.configure(source_dir, build_dir, effective_flags):
            logger.warning("Build configure failed—skipping auto-build (use pre-built)")
            return None  # Non-fatal: Fall back to command[0] as-is
        if not builder.build(build_dir, self.test.build_config.parallel_jobs):
            logger.warning("Build failed")
            return None
        
        # After successful build, auto-detect the actual binary name
        detected_binary = self._detect_binary_in_build_dir(build_dir)
        if detected_binary:
            # Update the executable name with what we actually found
            self.test.build_config.executable_name = detected_binary.name
            logger.info(f"Auto-detected binary: {detected_binary.name}")
        
        # Use the build directory - this is where the compiled binary actually is
        logger.info(f"Build completed - binary located in: {build_dir}")
        return build_dir
    
    def run(self) -> bool:
        logger.info(f"Starting {self.test.name} ({self.test.scaling_config.scaling_type.value})")
        
        job_configs = self.scaling_engine.generate_job_configs()
        logger.info(f"Generated {len(job_configs)} configs")
        
        # Debug: Log auto-submit status
        logger.info(f"Auto-submit enabled: {self.test.auto_submit_jobs}")
        
        job_ids = {}
        submission_failures = 0
        prepared_jobs = 0
        
        logger.info(f"\n{'='*60}")
        logger.info(f"AUTOMATED JOB SUBMISSION - {len(job_configs)} jobs to submit")
        logger.info(f"{'='*60}\n")
        
        for idx, job_config in enumerate(job_configs, 1):
            logger.info(f"[{idx}/{len(job_configs)}] Processing: {job_config.job_id}")
            logger.info(f"  Folder: output/.../{ job_config.job_id}/")
            
            try:
                job_id = self._submit_job(job_config)
                
                if job_id:
                    # Check if this is a prepared job (not actually submitted)
                    if job_id.startswith("prepared_"):
                        prepared_jobs += 1
                        logger.info(f"  Status: Prepared (not submitted)")
                    else:
                        job_ids[job_config.job_id] = job_id
                        logger.info(f"  Status: Submitted ✓")
                        logger.info(f"  Slurm ID: {job_id}\n")
                else:
                    submission_failures += 1
                    logger.error(f"  Status: Failed (returned None) ✗\n")
                    
            except subprocess.CalledProcessError as e:
                # sbatch command failed - this is CRITICAL
                submission_failures += 1
                logger.error(f"  Status: SUBMISSION FAILED ✗")
                logger.error(f"  Error: sbatch command failed with exit code {e.returncode}")
                logger.error(f"  stdout: {e.stdout}")
                logger.error(f"  stderr: {e.stderr}")
                # Don't continue if sbatch is broken
                if submission_failures >= 3:
                    logger.error(f"\n❌ STOPPING: Multiple submission failures detected")
                    logger.error(f"   This indicates sbatch is not working properly")
                    logger.error(f"   Please check:")
                    logger.error(f"   - Is Slurm installed? Run: which sbatch")
                    logger.error(f"   - Are you on the login node?")
                    logger.error(f"   - Check Slurm status: sinfo")
                    break
            except Exception as e:
                submission_failures += 1
                logger.error(f"  Status: Failed with exception ✗")
                logger.error(f"  Error: {e}")
                logger.error(f"  Type: {type(e).__name__}\n")
                # Stop after too many failures
                if submission_failures >= 3:
                    logger.error(f"\n❌ STOPPING: Too many submission failures")
                    break
        
        # Print submission summary
        logger.info(f"\n{'='*60}")
        logger.info(f"SUBMISSION SUMMARY")
        logger.info(f"{'='*60}")
        logger.info(f"  Total jobs: {len(job_configs)}")
        logger.info(f"  Submitted: {len(job_ids)}")
        logger.info(f"  Failed: {submission_failures}")
        if prepared_jobs > 0:
            logger.info(f"  Prepared (not submitted): {prepared_jobs}")
        logger.info(f"{'='*60}\n")
        
        # Handle case where auto-submit is disabled
        if not self.test.auto_submit_jobs:
            logger.info(f"All job scripts generated and saved. {prepared_jobs} jobs prepared for manual submission.")
            self.jobs_completed = False  # Jobs not submitted/completed
            return True
        
        # If all submissions failed, return False immediately
        if submission_failures == len(job_configs):
            logger.error("❌ All job submissions failed")
            return False
        elif submission_failures > 0:
            logger.warning(f"⚠ {submission_failures} out of {len(job_configs)} jobs failed to submit")
        
        # Only monitor jobs that were successfully submitted
        if job_ids:
            logger.info(f"📊 Monitoring {len(job_ids)} submitted jobs...")
            logger.info(f"   Polling interval: 10 seconds")
            logger.info(f"   Press Ctrl+C to stop monitoring and generate reports later\n")
            
            try:
                results = self._monitor_jobs(job_ids)
                self._generate_summary(job_configs, results)
                
                failures = [k for k, v in results.items() if v != JobStatus.COMPLETED]
                if failures:
                    logger.warning(f"⚠ Job failures: {failures}")
                    self.jobs_completed = False
                    return False
                    
                logger.info(f"\n{'='*60}")
                logger.info("✅ ALL JOBS COMPLETED SUCCESSFULLY")
                logger.info(f"{'='*60}\n")
                self.jobs_completed = True  # Mark jobs as completed
                return True
                
            except KeyboardInterrupt:
                logger.warning(f"\n⚠ Monitoring interrupted by user")
                logger.info(f"Jobs are still running in the background")
                logger.info(f"Check status with: squeue -u $USER")
                logger.info(f"Run report generation later when jobs complete")
                self.jobs_completed = False  # Jobs not completed yet
                return True  # Return success so workflow continues
                
        else:
            logger.error("No jobs were successfully submitted")
            return False

    def _submit_job(self, job_config: JobConfig) -> Optional[str]:
        import platform
        
        logger.debug(f"_submit_job called for {job_config.job_id}")
        
        # PRE-EXECUTION VALIDATION
        logger.info(f"\n{'='*60}")
        logger.info(f"VALIDATING: {job_config.job_id}")
        logger.info(f"{'='*60}")
        
        is_valid, errors = validate_job_config(job_config, self.test)
        if not is_valid:
            error_msg = f"\n❌ Configuration validation failed for {job_config.job_id}:\n"
            error_msg += "\n".join(f"  • {err}" for err in errors)
            logger.error(error_msg)
            raise ValidationError(error_msg)
        
        logger.info("✓ Validation passed")
        logger.info(f"  Decomposition: {job_config.procs_decomposition[0]}×{job_config.procs_decomposition[1]}×{job_config.procs_decomposition[2]} = {job_config.num_procs} procs")
        
        job_dir = (self.run_dir / job_config.job_id).resolve()  # Convert to absolute path
        create_directory(job_dir)
        job_config.working_dir = job_dir  # Now absolute path
        
        # Input file - make a copy of command to avoid modifying the original
        job_command = self.test.command.copy()
        
        if self.test.input_file:
            # Handle directory vs file input
            if self.test.input_file.is_dir():
                # Input_file is a directory - copy all files from it
                import shutil
                logger.info(f"Copying input directory: {self.test.input_file} -> {job_dir}/")
                
                # Find the main input file for validation
                input_files = sorted([f for f in self.test.input_file.iterdir() if f.is_file()])
                main_input = input_files[0] if input_files else None
                
                # Validate and adapt main input file if found
                if main_input:
                    logger.info(f"Validating and adapting input file: {main_input.name}")
                    adapted_input_path = job_dir / main_input.name
                    
                    # Perform validation and adaptation
                    base_procs = (self.test.scaling_config.initial_procs[0] * 
                                  self.test.scaling_config.initial_procs[1] * 
                                  self.test.scaling_config.initial_procs[2])
                    scale_factor = job_config.num_procs / base_procs if base_procs > 0 else 1.0
                    
                    try:
                        is_valid, output_file = validate_and_adapt_input(
                            input_file=main_input,
                            num_procs=job_config.num_procs,
                            procs_decomp=None,  # Auto-compute intelligent decomposition
                            scaling_type=self.test.scaling_config.scaling_type.value,
                            scale_factor=scale_factor,
                            output_file=adapted_input_path,
                            use_llm=self.test.use_llm_discovery,
                            llm_api_key=self.test.openai_api_key,
                            llm_model=self.test.llm_model,
                            custom_param_map=self.test.parameter_mapping
                        )
                        
                        if not is_valid:
                            logger.error(f"Input validation failed for {job_config.job_id}")
                            raise ValueError("Input file validation failed - see logs for details")
                        
                        logger.info(f"✓ Input file validated and adapted successfully")
                        
                    except Exception as e:
                        logger.warning(f"Input adaptation failed: {e}. Using original file.")
                        # Fall back to copying original
                        shutil.copy2(main_input, adapted_input_path)
                    
                    # Copy other files in directory
                    for item in input_files[1:]:
                        dest = job_dir / item.name
                        shutil.copy2(item, dest)
                        logger.debug(f"  Copied: {item.name}")
                    
                    input_name = main_input.name
                else:
                    # No files found, just copy directory
                    for item in self.test.input_file.iterdir():
                        if item.is_file():
                            dest = job_dir / item.name
                            shutil.copy2(item, dest)
                            logger.debug(f"  Copied: {item.name}")
                    input_name = input_files[0].name if input_files else "input.txt"
                
                # Update command with input filename
                if job_command:
                    last_arg = job_command[-1]
                    if isinstance(last_arg, str) and last_arg.lower().endswith((".inp", ".in", ".dat", ".txt", "stdin")):
                        job_command[-1] = input_name
                    else:
                        job_command.append(input_name)
                else:
                    job_command = [input_name]
            else:
                # Input_file is a single file - validate and adapt it
                original_name = self.test.input_file.name
                input_path = job_dir / original_name
                
                # Node 1 baseline: Use exact values from run.yaml to override placeholders
                if job_config.job_id == "node1":
                    logger.info("""╔═══════════════════════════════════════════════════════════════╗
║  NODE 1: Applying run.yaml values to input file           ║
╚═══════════════════════════════════════════════════════════════╝""")
                    
                    # Use generic input parser to apply run.yaml values
                    if self.test.scaling_config.variable_map:
                        logger.info("Using generic input parser for Node 1 (YAML-driven variable mapping)")
                        parser = GenericInputParser(self.test.scaling_config.variable_map)
                        
                        # Apply Node 1 values from run.yaml
                        success = parser.generate_scaled_input(
                            base_input=self.test.input_file,
                            output_file=input_path,
                            domain_size=self.test.scaling_config.initial_domain,
                            cell_count=self.test.scaling_config.initial_cells,
                            mpi_decomp=self.test.scaling_config.initial_procs,
                            particles_per_cell=self.test.scaling_config.particles_per_cell
                        )
                        
                        if success:
                            logger.info(f"✓ Node 1 input file generated with run.yaml values")
                        else:
                            logger.error("Failed to generate Node 1 input file, falling back to copy")
                            import shutil
                            shutil.copy2(self.test.input_file, input_path)
                    else:
                        # Fall back to copying if no variable map
                        import shutil
                        shutil.copy2(self.test.input_file, input_path)
                        logger.info(f"✓ Input file copied unchanged: {original_name}")
                
                else:
                    # For all other nodes: apply scaling/modifications
                    logger.info(f"Validating and adapting input file: {original_name}")
                    
                    # Calculate scale factor for weak scaling
                    base_procs = (self.test.scaling_config.initial_procs[0] * 
                                  self.test.scaling_config.initial_procs[1] * 
                                  self.test.scaling_config.initial_procs[2])
                    scale_factor = job_config.num_procs / base_procs if base_procs > 0 else 1.0
                    
                    # GENERIC PATH: Use GenericInputParser if variable_map is provided
                    if self.test.scaling_config.variable_map:
                        logger.info("Using generic input parser (YAML-driven variable mapping)")
                        parser = GenericInputParser(self.test.scaling_config.variable_map)
                        
                        try:
                            success = parser.generate_scaled_input(
                                base_input=self.test.input_file,
                                output_file=input_path,
                                domain_size=job_config.domain_size,
                                cell_count=job_config.cell_count,
                                mpi_decomp=job_config.procs_decomposition,
                                particles_per_cell=self.test.scaling_config.particles_per_cell
                            )
                            
                            if success:
                                logger.info(f"✓ Input file generated with generic parser")
                            else:
                                raise ValueError("Generic parser failed")
                        
                        except Exception as e:
                            logger.error(f"Generic parser failed: {e}")
                            # Fall back to copying original
                            import shutil
                            shutil.copy2(self.test.input_file, input_path)
                            logger.warning("Using unmodified input file as fallback")
                    
                    else:
                        # LLM-BASED PATH: Use existing LLM-based validation (if configured)
                        try:
                            is_valid, output_file = validate_and_adapt_input(
                                input_file=self.test.input_file,
                                num_procs=job_config.num_procs,
                                procs_decomp=None,  # Auto-compute intelligent decomposition
                                scaling_type=self.test.scaling_config.scaling_type.value,
                                scale_factor=scale_factor,
                                output_file=input_path,
                                use_llm=self.test.use_llm_discovery,
                                llm_api_key=self.test.openai_api_key,
                                llm_model=self.test.llm_model,
                                custom_param_map=self.test.parameter_mapping
                            )
                            
                            if not is_valid:
                                logger.error(f"Input validation failed for {job_config.job_id}")
                                raise ValueError("Input file validation failed - see logs for details")
                            
                            logger.info(f"✓ Input file validated and adapted successfully")
                            
                        except Exception as e:
                            logger.warning(f"Input adaptation failed: {e}. Using basic method.")
                            # Fall back to basic method
                            content = self.test.get_input_content(job_config)
                            write_file(input_path, content)
                
                # Update command with input filename
                if job_command:
                    last_arg = job_command[-1]
                    if isinstance(last_arg, str) and last_arg.lower().endswith((".inp", ".in", ".dat", ".txt", "stdin")):
                        job_command[-1] = original_name
                    else:
                        job_command.append(original_name)
                else:
                    # Should not happen due to validation, but be safe
                    job_command = [original_name]
        
        env_setup = self._generate_env_setup()
        launch_cmd = self.launcher.generate_launch_command(
            job_config, job_command, self.test.resource_config
        )
        
        # Apply QoS mapping based on node count for this specific job
        resource_config_for_job = self.test.resource_config
        selected_qos = resource_config_for_job.get_qos_for_nodes(job_config.num_nodes)
        if selected_qos and selected_qos != resource_config_for_job.qos:
            # Create a copy of resource config with the selected QoS
            import copy
            resource_config_for_job = copy.copy(resource_config_for_job)
            resource_config_for_job.qos = selected_qos
            logger.info(f"  QoS selected for {job_config.num_nodes} nodes: {selected_qos}")
        
        script_content = self.scheduler.generate_job_script(
            job_config, resource_config_for_job, launch_cmd, env_setup
        )
        
        # Use appropriate file extension based on platform and scheduler
        is_windows = platform.system() == "Windows"
        is_local_scheduler = self.test.backend_config.scheduler.name == "LOCAL"
        
        if is_windows and is_local_scheduler:
            script_path = job_dir / "job.bat"
        else:
            script_path = job_dir / "job.sh"
            
        script_path.write_text(script_content)
        if not (is_windows and is_local_scheduler):
            script_path.chmod(0o755)
        
        logger.info(f"Generated job script: {script_path}")
        logger.debug(f"Job script content:\n{script_content}")
        
        # Check if automatic job submission is enabled
        if not self.test.auto_submit_jobs:
            logger.info(f"Job script generated but not submitted (auto-submit disabled): {script_path}")
            # Create a metadata file to indicate the job was prepared but not submitted
            metadata = {
                'job_id': job_config.job_id,
                'num_nodes': job_config.num_nodes,
                'num_procs': job_config.num_procs,
                'procs_decomposition': job_config.procs_decomposition,
                'prepared_at': datetime.now().isoformat(),
                'script_path': str(script_path),
                'status': 'prepared'
            }
            with open(job_dir / "metadata.json", 'w') as f:
                json.dump(metadata, f, indent=2)
            return f"prepared_{job_config.job_id}"  # Return a special ID to indicate prepared but not submitted
        
        # Automatic submission is enabled, proceed with submission
        logger.info(f"▶ Submitting job to scheduler: {job_config.job_id}")
        logger.info(f"  Job script: {script_path}")
        logger.info(f"  Scheduler: {self.test.backend_config.scheduler.value}")
        logger.info(f"  Command: sbatch {script_path}")
        
        try:
            # This calls sbatch for Slurm scheduler
            job_id = self.scheduler.submit_job(script_path)
            
            if not job_id:
                raise ValueError("Scheduler returned empty job ID")
            
            metadata = {
                'job_id': job_id,
                'num_nodes': job_config.num_nodes,
                'num_procs': job_config.num_procs,
                'procs_decomposition': job_config.procs_decomposition,
                'submitted_at': datetime.now().isoformat(),
                'script_path': str(script_path),
                'status': 'submitted'
            }
            with open(job_dir / "metadata.json", 'w') as f:
                json.dump(metadata, f, indent=2)
            
            logger.info(f"✓ Successfully submitted job {job_id}")
            logger.info(f"  Slurm Job ID: {job_id}")
            logger.info(f"  Job folder: {job_dir}")
            return job_id
            
        except Exception as e:
            logger.error(f"✗ Submit failed for job {job_config.job_id}")
            logger.error(f"  Error: {e}")
            logger.error(f"  Script: {script_path}")
            logger.error(f"  Traceback:", exc_info=True)
            
            # Create metadata even for failed submissions
            metadata = {
                'job_id': job_config.job_id,
                'num_nodes': job_config.num_nodes,
                'num_procs': job_config.num_procs,
                'procs_decomposition': job_config.procs_decomposition,
                'prepared_at': datetime.now().isoformat(),
                'script_path': str(script_path),
                'status': 'submission_failed',
                'error': str(e)
            }
            with open(job_dir / "metadata.json", 'w') as f:
                json.dump(metadata, f, indent=2)
            
            # Re-raise to propagate the error
            raise
    
    def _generate_env_setup(self) -> List[str]:
        commands = []
        if self.test.environment_config.modules:
            commands.extend(self.module_system.generate_load_commands(self.test.environment_config.modules))
        for key, value in self.test.environment_config.env_vars.items():
            commands.append(f"export {key}={value}")
        commands.extend(self.test.environment_config.pre_commands)
        return commands
    
    def _monitor_jobs(self, job_ids: Dict[str, str]) -> Dict[str, JobStatus]:
        results = {}
        pending = set(job_ids.keys())
        while pending:
            for job_name in list(pending):
                job_id = job_ids[job_name]
                status = self.scheduler.get_job_status(job_id)
                if status in (JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED, JobStatus.TIMEOUT):
                    results[job_name] = status
                    pending.remove(job_name)
                    logger.info(f"{job_name}: {status.value}")
            if pending:
                time.sleep(10)
        return results
    
    def _generate_summary(self, job_configs: List[JobConfig], results: Dict[str, JobStatus]):
        summary = {
            'test_name': self.test.name,
            'scaling_type': self.test.scaling_config.scaling_type.value,
            'completed_at': datetime.now().isoformat(),
            'jobs': []
        }
        parser = DefaultResultParser()
        for config in job_configs:
            job_dir = self.run_dir / config.job_id
            job_info = {
                'job_id': config.job_id,
                'num_nodes': config.num_nodes,
                'num_procs': config.num_procs,
                'procs_decomposition': list(config.procs_decomposition),  # Add decomposition
                'status': results.get(config.job_id, JobStatus.UNKNOWN).value
            }
            
            # Try multiple output file patterns (in priority order)
            out_files = (
                list(job_dir.glob("job.out")) or
                list(job_dir.glob("slurm-*.out")) or
                list(job_dir.glob("out_*.out")) or
                list(job_dir.glob("*.out"))
            )
            
            # Also check for timing.json
            timing_json = job_dir / "timing.json"
            if timing_json.exists():
                try:
                    with open(timing_json, 'r') as f:
                        timing_data = json.load(f)
                    if 'wall_time' in timing_data:
                        job_info['wall_time'] = timing_data['wall_time']
                        job_info['metrics'] = timing_data
                        logger.info(f"  Loaded timing from {timing_json.name}: {timing_data.get('wall_time')}s")
                except Exception as e:
                    logger.warning(f"  Failed to parse timing.json in {job_dir.name}: {e}")
            
            # Fallback: parse from output files
            if 'wall_time' not in job_info and out_files:
                logger.debug(f"  Parsing timing from output file: {out_files[0].name}")
                try:
                    metrics = parser.parse_output(out_files[0])
                    job_info['metrics'] = metrics
                    if 'wall_time' in metrics:
                        job_info['wall_time'] = metrics['wall_time']
                        logger.info(f"  Extracted timing from {out_files[0].name}: {metrics['wall_time']}s")
                except Exception as e:
                    logger.warning(f"  Failed to parse output file {out_files[0].name}: {e}")
            
            if 'wall_time' not in job_info:
                logger.warning(f"  No timing data found for {config.job_id} in {job_dir}")
            
            summary['jobs'].append(job_info)
        
        summary_path = self.run_dir / "summary.json"
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        logger.info(f"Summary: {summary_path}")
