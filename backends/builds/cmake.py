"""
Enhanced CMake backend with NVHPC support and CUDA library path handling - BUGFIXED

CRITICAL FIXES:
- Added NVHPC detection and special handling
- Added CUDA library path configuration
- Added CMAKE_CUDA_COMPILER_FORCED for NVHPC compatibility
- Added CMAKE_CUDA_HOST_COMPILER configuration
- Improved error messages for CUDA-related failures

This module:
1. ACTUALLY runs cmake and make commands (not just logs them)
2. Handles NVHPC's bundled CUDA toolkit properly
3. Sets correct CUDA library paths
4. Verifies that compilation produces output files
5. Automatically detects the compiled binary name
"""

import subprocess
import logging
import shutil
import stat
import os
import sys
import shlex
from pathlib import Path
from typing import Optional, Dict, List, Tuple

from core.abstracts import BuildSystemInterface


logger = logging.getLogger(__name__)


class CMakeBackend(BuildSystemInterface):
    """CMake backend with NVHPC support and proper CUDA handling."""
    
    def __init__(self, options: Optional[Dict] = None):
        """Initialize CMake backend.
        
        Args:
            options: Dict with optional 'module_commands' list
        """
        self.options = options or {}
        self.module_commands = self.options.get('module_commands', [])
        self.detected_executable: Optional[Path] = None
        self._nvhpc_detected = False
        self._cuda_root = None
        
        # Detect NVHPC from module commands
        self._detect_nvhpc()
        
        # Log initialization
        logger.info(f"CMakeBackend initialized")
        if self.module_commands:
            logger.info(f"  Module commands: {len(self.module_commands)}")
            for cmd in self.module_commands:
                logger.info(f"    - {cmd}")
        if self._nvhpc_detected:
            logger.info(f"  NVHPC detected: Will apply NVHPC-specific configurations")
    
    def _detect_nvhpc(self):
        """Detect if NVHPC is being used from module commands."""
        for cmd in self.module_commands:
            if 'nvhpc' in cmd.lower():
                self._nvhpc_detected = True
                logger.debug("NVHPC module detected")
                break
    
    def _find_cuda_root(self) -> Optional[str]:
        """
        Find CUDA installation root with improved detection.
        
        Priority:
        1. Module-loaded CUDA (CUDA_HOME, CUDA_ROOT from loaded modules)
        2. NVHPC bundled CUDA
        3. Common system locations
        
        Returns:
            Path to CUDA root or None
        """
        if self._cuda_root:
            return self._cuda_root
        
        # Priority 1: Check environment variables set by loaded modules
        for var in ['CUDA_HOME', 'CUDA_ROOT', 'CUDA_PATH', 'CUDADIR']:
            if var in os.environ:
                cuda_path = os.environ[var]
                if os.path.exists(cuda_path):
                    # Verify this is a valid CUDA installation
                    if self._validate_cuda_path(cuda_path):
                        self._cuda_root = cuda_path
                        logger.info(f"Found CUDA root from ${var}: {cuda_path}")
                        return self._cuda_root
        
        # Priority 2: For NVHPC, try to find bundled CUDA
        if self._nvhpc_detected:
            # Check NVHPC environment variables
            for nvhpc_var in ['NVHPC_ROOT', 'NVHPC_HOME', 'NVHPC_DIR']:
                nvhpc_root = os.environ.get(nvhpc_var, '')
                if nvhpc_root:
                    # NVHPC bundles CUDA in Linux_x86_64/<version>/cuda/
                    cuda_candidates = [
                        Path(nvhpc_root) / "Linux_x86_64" / "cuda",
                        Path(nvhpc_root) / "cuda",
                    ]
                    for candidate in cuda_candidates:
                        if candidate.exists() and self._validate_cuda_path(str(candidate)):
                            self._cuda_root = str(candidate)
                            logger.info(f"Found NVHPC bundled CUDA: {self._cuda_root}")
                            return self._cuda_root
        
        # Priority 3: Try common system locations
        common_paths = [
            '/usr/local/cuda',
            '/opt/nvidia/hpc_sdk/Linux_x86_64/cuda',  # NVHPC SDK location
            '/opt/cuda',
            '/usr/cuda',
            '/usr/local/cuda-12.6',  # Version-specific paths
            '/usr/local/cuda-12.4',
            '/usr/local/cuda-12',
        ]
        for path in common_paths:
            if os.path.exists(path) and self._validate_cuda_path(path):
                self._cuda_root = path
                logger.info(f"Found CUDA at common location: {path}")
                return self._cuda_root
        
        logger.warning("Could not find valid CUDA root directory")
        return None
    
    def _validate_cuda_path(self, cuda_path: str) -> bool:
        """
        Validate that a path contains a valid CUDA installation.
        
        Args:
            cuda_path: Path to check
            
        Returns:
            True if valid CUDA installation
        """
        # Check for essential CUDA components
        nvcc = os.path.join(cuda_path, 'bin', 'nvcc')
        
        # Check for library directories
        has_libs = False
        for lib_dir in ['lib64', 'lib', 'lib/x86_64-linux-gnu']:
            lib_path = os.path.join(cuda_path, lib_dir)
            if os.path.exists(lib_path):
                # Check for essential CUDA libraries
                cudart = os.path.join(lib_path, 'libcudart.so')
                cudart_static = os.path.join(lib_path, 'libcudart_static.a')
                if os.path.exists(cudart) or os.path.exists(cudart_static):
                    has_libs = True
                    break
        
        return os.path.exists(nvcc) and has_libs
    
    def _add_nvhpc_cuda_flags(self, flags: Dict[str, str]) -> Dict[str, str]:
        """
        Add NVHPC-specific CUDA flags.
        
        Args:
            flags: Existing CMake flags
            
        Returns:
            Updated flags dict
        """
        if not self._nvhpc_detected:
            return flags
        
        logger.info("Adding NVHPC-specific CUDA configurations...")
        
        # Force CMake to accept the CUDA compiler without testing
        # This avoids the compute_52 default test that fails with NVHPC
        if 'CMAKE_CUDA_COMPILER_FORCED' not in flags:
            flags['CMAKE_CUDA_COMPILER_FORCED'] = 'ON'
            logger.info("  Set CMAKE_CUDA_COMPILER_FORCED=ON")
        
        # Set CUDA host compiler explicitly
        if 'CMAKE_CUDA_HOST_COMPILER' not in flags:
            # Use the C++ compiler as the CUDA host compiler
            cxx_compiler = os.environ.get('CXX', 'g++')
            flags['CMAKE_CUDA_HOST_COMPILER'] = cxx_compiler
            logger.info(f"  Set CMAKE_CUDA_HOST_COMPILER={cxx_compiler}")
        
        # Find and set CUDA root
        cuda_root = self._find_cuda_root()
        if cuda_root:
            if 'CUDA_TOOLKIT_ROOT_DIR' not in flags:
                flags['CUDA_TOOLKIT_ROOT_DIR'] = cuda_root
                logger.info(f"  Set CUDA_TOOLKIT_ROOT_DIR={cuda_root}")
            
            if 'CMAKE_CUDA_COMPILER' not in flags:
                nvcc_path = os.path.join(cuda_root, 'bin', 'nvcc')
                if os.path.exists(nvcc_path):
                    flags['CMAKE_CUDA_COMPILER'] = nvcc_path
                    logger.info(f"  Set CMAKE_CUDA_COMPILER={nvcc_path}")
            
            # Determine library directory paths (try multiple common locations)
            cuda_lib_dirs = []
            for lib_subdir in ['lib64', 'lib', 'lib/x86_64-linux-gnu']:
                lib_path = os.path.join(cuda_root, lib_subdir)
                if os.path.exists(lib_path):
                    cuda_lib_dirs.append(lib_path)
            
            if cuda_lib_dirs:
                # Add primary library directory to CMAKE_PREFIX_PATH
                existing_prefix = flags.get('CMAKE_PREFIX_PATH', '')
                if existing_prefix:
                    flags['CMAKE_PREFIX_PATH'] = f"{cuda_root};{existing_prefix}"
                else:
                    flags['CMAKE_PREFIX_PATH'] = cuda_root
                logger.info(f"  Set CMAKE_PREFIX_PATH to include CUDA root")
                
                # CRITICAL FIX: Explicitly set CUDA library directories
                # This ensures the linker can find libcudadevrt.a and libcudart_static.a
                cuda_lib_dirs_str = ';'.join(cuda_lib_dirs)
                flags['CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES'] = cuda_lib_dirs_str
                logger.info(f"  Set CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES={cuda_lib_dirs_str}")
                
                # Also set CUDA runtime library paths explicitly
                flags['CUDA_CUDART_LIBRARY'] = os.path.join(cuda_lib_dirs[0], 'libcudart.so')
                logger.info(f"  Set CUDA_CUDART_LIBRARY")
                
                # Set linker flags to include CUDA library paths
                # This is crucial for NVHPC + CUDA linking
                link_flags = []
                for lib_dir in cuda_lib_dirs:
                    link_flags.append(f"-L{lib_dir}")
                link_flags.append("-lcudart")
                link_flags.append("-lcudadevrt")
                
                cuda_link_flags = ' '.join(link_flags)
                
                # Add to CMAKE_EXE_LINKER_FLAGS
                existing_linker_flags = flags.get('CMAKE_EXE_LINKER_FLAGS', '')
                if existing_linker_flags:
                    flags['CMAKE_EXE_LINKER_FLAGS'] = f"{existing_linker_flags} {cuda_link_flags}"
                else:
                    flags['CMAKE_EXE_LINKER_FLAGS'] = cuda_link_flags
                
                # Add to CMAKE_SHARED_LINKER_FLAGS (for shared libraries)
                existing_shared_flags = flags.get('CMAKE_SHARED_LINKER_FLAGS', '')
                if existing_shared_flags:
                    flags['CMAKE_SHARED_LINKER_FLAGS'] = f"{existing_shared_flags} {cuda_link_flags}"
                else:
                    flags['CMAKE_SHARED_LINKER_FLAGS'] = cuda_link_flags
                
                logger.info(f"  Set CMAKE_EXE_LINKER_FLAGS with CUDA library paths")
                logger.info(f"  Set CMAKE_SHARED_LINKER_FLAGS with CUDA library paths")
            else:
                logger.warning(f"  No CUDA library directories found in {cuda_root}")
        
        return flags
    
    def _execute_command(
        self, 
        cmd: List[str], 
        cwd: Path, 
        description: str,
        timeout: int = 3600
    ) -> Tuple[bool, str, str]:
        """
        Execute a command and return results.
        
        This method ACTUALLY runs the command (not just logs it).
        
        Args:
            cmd: Command as list of strings
            cwd: Working directory
            description: Human-readable description
            timeout: Timeout in seconds
            
        Returns:
            Tuple of (success, stdout, stderr)
        """
        cmd_str = ' '.join(cmd)
        
        logger.info(f"")
        logger.info(f">>> EXECUTING: {description}")
        logger.info(f">>> Directory: {cwd}")
        
        try:
            if self.module_commands:
                # Build full command with module loads
                full_cmd = " && ".join(self.module_commands + [cmd_str])
                shell_cmd = f"bash -l -c '{full_cmd}'"
                logger.info(f">>> Command (with modules): {shell_cmd[:200]}...")
                
                result = subprocess.run(
                    shell_cmd,
                    shell=True,
                    cwd=str(cwd),
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                    env=os.environ.copy()
                )
            else:
                logger.info(f">>> Command: {cmd_str}")
                
                result = subprocess.run(
                    cmd,
                    cwd=str(cwd),
                    capture_output=True,
                    text=True,
                    timeout=timeout,
                    env=os.environ.copy()
                )
            
            success = (result.returncode == 0)
            
            if success:
                logger.info(f">>> SUCCESS (exit code 0)")
                if result.stdout:
                    # Log last few lines of output
                    lines = result.stdout.strip().split('\n')
                    if len(lines) > 5:
                        logger.info(f">>> Output (last 5 lines):")
                        for line in lines[-5:]:
                            logger.info(f"    {line}")
                    else:
                        for line in lines:
                            logger.info(f"    {line}")
            else:
                logger.error(f">>> FAILED (exit code {result.returncode})")
                
                # Enhanced CUDA-specific error detection
                stderr_lower = result.stderr.lower()
                stdout_lower = result.stdout.lower()
                combined = stderr_lower + stdout_lower
                
                cuda_error_detected = False
                
                if any(term in combined for term in ['cuda', 'nvcc', 'cudart', 'cudadevrt']):
                    cuda_error_detected = True
                    logger.error(f"")
                    logger.error(f"CUDA-RELATED ERROR DETECTED:")
                    
                    # Check for specific CUDA library linking errors
                    if 'cannot find -lcudart' in result.stderr or 'cannot find -lcudadevrt' in result.stderr:
                        logger.error(f"  - CUDA runtime libraries not found during linking")
                        logger.error(f"  - Missing: libcudart_static.a and/or libcudadevrt.a")
                        logger.error(f"  - This typically means CUDA library paths are not in linker search path")
                        logger.error(f"")
                        logger.error(f"  RESOLUTION STEPS:")
                        logger.error(f"    1. Verify CUDA is properly installed")
                        logger.error(f"    2. Check that CUDA libraries exist in:")
                        cuda_root = self._cuda_root or self._find_cuda_root()
                        if cuda_root:
                            for lib_dir in ['lib64', 'lib']:
                                check_path = os.path.join(cuda_root, lib_dir)
                                if os.path.exists(check_path):
                                    logger.error(f"       {check_path}")
                        logger.error(f"    3. Ensure CMAKE_SHARED_LINKER_FLAGS includes CUDA library paths")
                        logger.error(f"    4. If using NVHPC, ensure CUDA module is loaded")
                    
                    # Check for CUDA architecture errors
                    if 'compute_' in result.stderr:
                        logger.error(f"  - Wrong CUDA architecture specified")
                        logger.error(f"  - Check CMAKE_CUDA_ARCHITECTURES setting")
                    
                    # Check for general CUDA library path issues
                    if 'cannot find -lcuda' in result.stderr and 'cudart' not in result.stderr:
                        logger.error(f"  - CUDA driver library not found")
                        logger.error(f"  - Check LD_LIBRARY_PATH and CUDA installation")
                    
                    # Check for CUDA compiler issues
                    if 'nvcc' in combined and 'not found' in combined:
                        logger.error(f"  - CUDA compiler (nvcc) not found")
                        logger.error(f"  - Verify CUDA_TOOLKIT_ROOT_DIR is set correctly")
                    
                    logger.error(f"")
                
                if result.stderr:
                    logger.error(f">>> STDERR:")
                    for line in result.stderr.strip().split('\n')[-20:]:
                        logger.error(f"    {line}")
                if result.stdout:
                    logger.error(f">>> STDOUT:")
                    for line in result.stdout.strip().split('\n')[-10:]:
                        logger.error(f"    {line}")
            
            return success, result.stdout, result.stderr
            
        except subprocess.TimeoutExpired:
            logger.error(f">>> TIMEOUT after {timeout} seconds")
            return False, "", f"Command timed out after {timeout}s"
        except FileNotFoundError as e:
            logger.error(f">>> COMMAND NOT FOUND: {e}")
            return False, "", str(e)
        except Exception as e:
            logger.error(f">>> EXCEPTION: {type(e).__name__}: {e}")
            return False, "", str(e)
    
    def _scan_for_binaries(self, directory: Path) -> List[Path]:
        """
        Scan directory for compiled binary executables.
        
        Identifies ELF binaries by reading file headers.
        No filename assumptions - purely based on file content.
        
        Args:
            directory: Directory to scan
            
        Returns:
            List of paths to binary executables
        """
        binaries = []
        
        # Files/directories to skip
        skip_names = {
            'CMakeFiles', 'cmake_install.cmake', 'CMakeCache.txt',
            'Makefile', 'CTestTestfile.cmake', 'compile_commands.json',
            'CPackConfig.cmake', 'CPackSourceConfig.cmake',
            '.ninja_deps', '.ninja_log', 'build.ninja', 'rules.ninja'
        }
        
        skip_extensions = {
            '.o', '.a', '.so', '.la', '.lo', '.dylib', '.dll',
            '.cmake', '.txt', '.md', '.h', '.hpp', '.c', '.cpp',
            '.cu', '.f', '.f90', '.py', '.sh', '.pl', '.rb',
            '.json', '.xml', '.yaml', '.yml', '.in', '.inc'
        }
        
        def scan_dir(path: Path, depth: int = 0):
            if depth > 3:  # Limit recursion
                return
            
            try:
                for item in path.iterdir():
                    if item.name in skip_names or item.name.startswith('.'):
                        continue
                    
                    if item.is_dir():
                        if item.name not in {'CMakeFiles', '__pycache__', '.git'}:
                            scan_dir(item, depth + 1)
                    
                    elif item.is_file():
                        # Skip by extension
                        if item.suffix.lower() in skip_extensions:
                            continue
                        
                        # Skip library-like names
                        if item.name.startswith('lib'):
                            continue
                        
                        # Check if executable
                        try:
                            st = item.stat()
                            if not (st.st_mode & stat.S_IXUSR):
                                continue
                            
                            # Read file header to verify it's a binary
                            with open(item, 'rb') as f:
                                header = f.read(4)
                            
                            # ELF magic number
                            if header[:4] == b'\x7fELF':
                                binaries.append(item)
                                logger.debug(f"Found ELF binary: {item}")
                            # Mach-O (macOS)
                            elif header[:4] in (b'\xfe\xed\xfa\xce', b'\xfe\xed\xfa\xcf',
                                               b'\xca\xfe\xba\xbe', b'\xcf\xfa\xed\xfe'):
                                binaries.append(item)
                                logger.debug(f"Found Mach-O binary: {item}")
                                
                        except (OSError, IOError):
                            pass
                            
            except PermissionError:
                pass
        
        scan_dir(directory)
        return binaries
    
    def _choose_best_binary(
        self, 
        binaries: List[Path], 
        source_dir: Optional[Path] = None
    ) -> Optional[Path]:
        """
        Choose the most likely main executable from candidates.
        
        Uses heuristics based on name matching and location.
        
        Args:
            binaries: List of candidate binaries
            source_dir: Source directory for name hints
            
        Returns:
            Best candidate or None
        """
        if not binaries:
            return None
        
        if len(binaries) == 1:
            return binaries[0]
        
        # Get hint from source directory name
        hint = ""
        if source_dir:
            hint = source_dir.name.lower()
            # Remove common suffixes to get core name
            for suffix in ['-cpu', '-gpu', '-mpi', '-omp', '-ns', '-cuda',
                          '_cpu', '_gpu', '_mpi', '_omp', '_ns', '_cuda']:
                hint = hint.replace(suffix, '')
        
        # Score each binary
        scored = []
        for binary in binaries:
            score = 0
            name = binary.name.lower()
            stem = binary.stem.lower()
            
            # Exact match with hint
            if hint and stem == hint:
                score += 50
            
            # Hint contained in name
            if hint and hint in name:
                score += 30
            
            # Name contained in hint (e.g., "ipic3d" in "ipic3d-cpu-ns")
            if hint:
                clean_hint = hint.replace('-', '').replace('_', '')
                clean_name = stem.replace('-', '').replace('_', '')
                if clean_name in clean_hint or clean_hint in clean_name:
                    score += 20
            
            # Penalize test/example executables
            if any(word in name for word in ['test', 'example', 'demo', 'bench']):
                score -= 30
            
            # Prefer bin/ directory
            if 'bin' in str(binary.parent):
                score += 15
            
            # Prefer shorter names (main executables tend to have simple names)
            score -= len(name) * 0.1
            
            scored.append((score, binary))
        
        # Sort by score (descending)
        scored.sort(key=lambda x: x[0], reverse=True)
        
        logger.debug(f"Binary scores:")
        for score, binary in scored[:5]:
            logger.debug(f"  {score:6.1f} - {binary.name}")
        
        return scored[0][1]
    
    def configure(
        self, 
        source_dir: Path, 
        build_dir: Path, 
        flags: Optional[Dict[str, str]] = None
    ) -> bool:
        """
        Run cmake configuration.
        
        Args:
            source_dir: Path to source code
            build_dir: Path to build directory
            flags: CMake -D flags
            
        Returns:
            True if cmake succeeded
        """
        source_dir = Path(source_dir).resolve()
        build_dir = Path(build_dir).resolve()
        
        # Create build directory
        build_dir.mkdir(parents=True, exist_ok=True)
        
        # Apply NVHPC-specific flags if needed
        effective_flags = dict(flags or {})
        if self._nvhpc_detected:
            effective_flags = self._add_nvhpc_cuda_flags(effective_flags)
        
        # Build cmake command
        cmd = ["cmake", str(source_dir)]
        if effective_flags:
            for key, value in effective_flags.items():
                # Quote values that contain spaces to prevent shell from splitting them
                # Don't quote semicolons - they're CMake's list separator and should stay unquoted
                if ' ' in str(value):
                    cmd.append(f'-D{key}="{value}"')
                else:
                    cmd.append(f"-D{key}={value}")
        
        logger.info(f"")
        logger.info(f"{'='*60}")
        logger.info(f"STEP 1: CMAKE CONFIGURE")
        logger.info(f"{'='*60}")
        logger.info(f"Source:  {source_dir}")
        logger.info(f"Build:   {build_dir}")
        logger.info(f"Flags:   {effective_flags}")
        
        success, stdout, stderr = self._execute_command(
            cmd, build_dir, "cmake configure"
        )
        
        if success:
            # Verify CMakeCache.txt was created
            cache = build_dir / "CMakeCache.txt"
            if cache.exists():
                logger.info(f"✓ CMakeCache.txt created")
            else:
                logger.warning(f"⚠ CMakeCache.txt not found - configure may have failed")
                success = False
        
        return success
    
    def build(self, build_dir: Path, parallel_jobs: int = 4) -> bool:
        """
        Run the actual build (make).
        
        Args:
            build_dir: Path to build directory
            parallel_jobs: Number of parallel jobs
            
        Returns:
            True if build succeeded
        """
        build_dir = Path(build_dir).resolve()
        
        logger.info(f"")
        logger.info(f"{'='*60}")
        logger.info(f"STEP 2: BUILD (COMPILE)")
        logger.info(f"{'='*60}")
        logger.info(f"Build dir:      {build_dir}")
        logger.info(f"Parallel jobs:  {parallel_jobs}")
        
        # Try cmake --build first
        cmd = ["cmake", "--build", ".", "--parallel", str(parallel_jobs)]
        success, stdout, stderr = self._execute_command(
            cmd, build_dir, "cmake --build"
        )
        
        if not success:
            # Fallback to make
            logger.info(f"Trying fallback: make -j{parallel_jobs}")
            cmd = ["make", f"-j{parallel_jobs}"]
            success, stdout, stderr = self._execute_command(
                cmd, build_dir, "make"
            )
        
        if success:
            logger.info(f"✓ Build completed")
        else:
            logger.error(f"✗ Build FAILED")
        
        return success
    
    def find_executable(
        self, 
        build_dir: Path, 
        source_dir: Optional[Path] = None
    ) -> Optional[Path]:
        """
        Find the compiled executable in build directory.
        
        Scans for ELF binaries and selects the best match.
        NO hardcoded names - purely detection based.
        
        Args:
            build_dir: Path to build directory
            source_dir: Path to source (for name hints)
            
        Returns:
            Path to executable or None
        """
        build_dir = Path(build_dir).resolve()
        
        logger.info(f"")
        logger.info(f"{'='*60}")
        logger.info(f"STEP 3: DETECT COMPILED BINARY")
        logger.info(f"{'='*60}")
        logger.info(f"Scanning: {build_dir}")
        
        # Scan for binaries
        binaries = self._scan_for_binaries(build_dir)
        
        if not binaries:
            logger.error(f"✗ NO BINARIES FOUND in build directory!")
            logger.error(f"  Build directory contents:")
            try:
                for item in sorted(build_dir.iterdir())[:20]:
                    item_type = "[DIR]" if item.is_dir() else "[FILE]"
                    logger.error(f"    {item_type} {item.name}")
            except Exception as e:
                logger.error(f"    Error listing: {e}")
            return None
        
        logger.info(f"Found {len(binaries)} binary file(s)")
        
        # Choose the best one
        self.detected_executable = self._choose_best_binary(binaries, source_dir)
        
        if self.detected_executable:
            logger.info(f"")
            logger.info(f"✓ DETECTED EXECUTABLE: {self.detected_executable.name}")
            logger.info(f"✓ Full path: {self.detected_executable}")
        
        return self.detected_executable
    
    def build_and_find(
        self, 
        source_dir: Path, 
        build_dir: Path,
        flags: Optional[Dict[str, str]] = None,
        parallel_jobs: int = 4,
        force_rebuild: bool = False
    ) -> Optional[Path]:
        """
        Complete build workflow: configure, compile, detect binary.
        
        This is the MAIN ENTRY POINT.
        
        After calling this:
        - self.detected_executable contains the binary path
        - Returns the path to the compiled executable
        
        Args:
            source_dir: Path to source code
            build_dir: Path for build output
            flags: CMake flags
            parallel_jobs: Parallel build jobs
            force_rebuild: If True, rebuild even if binary exists
            
        Returns:
            Path to compiled executable, or None if failed
        """
        source_dir = Path(source_dir).resolve()
        build_dir = Path(build_dir).resolve()
        
        logger.info(f"")
        logger.info(f"{'#'*60}")
        logger.info(f"#  AUTOMATIC BUILD AND BINARY DETECTION")
        logger.info(f"{'#'*60}")
        logger.info(f"#")
        logger.info(f"#  Source:    {source_dir}")
        logger.info(f"#  Build:     {build_dir}")
        logger.info(f"#  Jobs:      {parallel_jobs}")
        logger.info(f"#")
        logger.info(f"{'#'*60}")
        logger.info(f"")
        
        # BINARY CACHING: Check if build already exists and is recent
        if not force_rebuild and build_dir.exists():
            logger.info(f"✓ Build directory exists, checking for existing binary...")
            
            # Try to find existing executable
            existing_binary = self.find_executable(build_dir, source_dir)
            
            if existing_binary and existing_binary.exists():
                # Check if binary is newer than source files
                binary_time = existing_binary.stat().st_mtime
                
                # Check a few key source files
                source_newer = False
                cmake_file = source_dir / "CMakeLists.txt"
                if cmake_file.exists() and cmake_file.stat().st_mtime > binary_time:
                    source_newer = True
                
                if not source_newer:
                    logger.info(f"")
                    logger.info(f"{'#'*60}")
                    logger.info(f"#  USING CACHED BINARY (skipping rebuild)")
                    logger.info(f"{'#'*60}")
                    logger.info(f"#")
                    logger.info(f"#  Binary:    {existing_binary.name}")
                    logger.info(f"#  Location:  {existing_binary.parent}")
                    logger.info(f"#  Built:     {existing_binary.stat().st_mtime}")
                    logger.info(f"#")
                    logger.info(f"#  To force rebuild, delete build directory or use --rebuild flag")
                    logger.info(f"#")
                    logger.info(f"{'#'*60}")
                    logger.info(f"")
                    
                    self.detected_executable = existing_binary
                    return existing_binary
                else:
                    logger.info(f"⚠ Source files newer than binary, rebuilding...")
            else:
                logger.info(f"⚠ No existing binary found, building...")
        
        # Proceed with normal build
        logger.info(f"")
        logger.info(f"Starting build process...")
        logger.info(f"")
        
        # Step 1: Configure
        if not self.configure(source_dir, build_dir, flags):
            logger.error(f"")
            logger.error(f"{'!'*60}")
            logger.error(f"BUILD FAILED AT CONFIGURE STEP")
            logger.error(f"{'!'*60}")
            return None
        
        # Step 2: Build
        if not self.build(build_dir, parallel_jobs):
            logger.error(f"")
            logger.error(f"{'!'*60}")
            logger.error(f"BUILD FAILED AT COMPILE STEP")
            logger.error(f"{'!'*60}")
            return None
        
        # Step 3: Find executable
        executable = self.find_executable(build_dir, source_dir)
        
        if executable:
            logger.info(f"")
            logger.info(f"{'#'*60}")
            logger.info(f"#  BUILD SUCCESSFUL")
            logger.info(f"{'#'*60}")
            logger.info(f"#")
            logger.info(f"#  Executable: {executable.name}")
            logger.info(f"#  Location:   {executable.parent}")
            logger.info(f"#")
            logger.info(f"{'#'*60}")
            logger.info(f"")
        else:
            logger.error(f"")
            logger.error(f"{'!'*60}")
            logger.error(f"BUILD COMPLETED BUT NO BINARY FOUND")
            logger.error(f"{'!'*60}")
        
        return executable
    
    def install(self, build_dir: Path, install_dir: Path) -> bool:
        """Run cmake install."""
        build_dir = Path(build_dir).resolve()
        install_dir = Path(install_dir).resolve()
        
        cmd = ["cmake", "--install", ".", "--prefix", str(install_dir)]
        success, _, _ = self._execute_command(cmd, build_dir, "cmake install")
        return success
    
    def clean(self, build_dir: Path) -> bool:
        """Clean the build directory."""
        build_dir = Path(build_dir).resolve()
        
        cmd = ["cmake", "--build", ".", "--target", "clean"]
        success, _, _ = self._execute_command(cmd, build_dir, "cmake clean")
        
        if not success:
            # Manual cleanup
            logger.info("Manual cleanup of build directory")
            try:
                for item in build_dir.iterdir():
                    if item.is_dir():
                        shutil.rmtree(item)
                    else:
                        item.unlink()
                success = True
            except Exception as e:
                logger.error(f"Cleanup failed: {e}")
        
        return success
