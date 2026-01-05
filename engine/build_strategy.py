"""
Build Strategy Module - Automated Build Configuration Selection

This module is responsible for:
1. Detecting application requirements (CPU/GPU/hybrid)
2. Detecting available system hardware
3. Selecting appropriate build strategy
4. Configuring CMake/build flags accordingly

Design Principles:
- FULLY AUTOMATED: No manual intervention required
- SYSTEM AGNOSTIC: Works on any HPC system
- COMPILER AGNOSTIC: Detects available compilers
- APPLICATION AGNOSTIC: Analyzes source code to determine requirements

Author: HPC-ScaleTest Contributors
"""

import logging
import subprocess
import os
import re
import shutil
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple
from enum import Enum

logger = logging.getLogger(__name__)


class ExecutionModel(Enum):
    """Execution model for the application."""
    CPU_ONLY = "cpu"
    GPU_ONLY = "gpu"
    HYBRID = "hybrid"  # CPU + GPU acceleration
    UNKNOWN = "unknown"


class BuildType(Enum):
    """Build configuration type."""
    CPU_RELEASE = "cpu_release"
    GPU_CUDA = "gpu_cuda"
    GPU_HIP = "gpu_hip"        # AMD ROCm
    GPU_SYCL = "gpu_sycl"      # Intel oneAPI
    HYBRID_CUDA = "hybrid_cuda"
    HYBRID_HIP = "hybrid_hip"


@dataclass
class ApplicationRequirements:
    """
    Detected application requirements from source code analysis.
    
    All fields are detected at runtime - no hardcoded assumptions.
    """
    # Parallel programming models
    requires_mpi: bool = False
    requires_openmp: bool = False
    requires_cuda: bool = False
    requires_hip: bool = False      # AMD ROCm
    requires_sycl: bool = False     # Intel oneAPI
    requires_openacc: bool = False
    
    # Libraries
    requires_hdf5: bool = False
    requires_fftw: bool = False
    requires_petsc: bool = False
    requires_blas_lapack: bool = False
    
    # Build system
    build_system: str = "cmake"  # cmake, make, autotools
    
    # Detected files
    cuda_files: List[str] = field(default_factory=list)
    hip_files: List[str] = field(default_factory=list)
    source_files: List[str] = field(default_factory=list)
    
    # Recommended execution model
    execution_model: ExecutionModel = ExecutionModel.UNKNOWN
    
    # Confidence score (0-1)
    confidence: float = 0.0


@dataclass
class SystemCapabilities:
    """
    Detected system hardware capabilities.
    
    All fields are detected at runtime from SLURM or system queries.
    """
    # CPU
    cpu_cores_per_node: Optional[int] = None
    cpu_model: Optional[str] = None
    
    # GPU
    gpus_per_node: int = 0
    gpu_vendor: Optional[str] = None  # nvidia, amd, intel
    gpu_model: Optional[str] = None
    gpu_architecture: Optional[str] = None  # sm_80, gfx90a, etc.
    
    # Memory
    memory_per_node_gb: Optional[float] = None
    gpu_memory_gb: Optional[float] = None
    
    # Available compilers (detected from modules or PATH)
    available_compilers: List[str] = field(default_factory=list)
    available_cuda_compilers: List[str] = field(default_factory=list)
    available_hip_compilers: List[str] = field(default_factory=list)
    
    # Partition info
    partition: Optional[str] = None
    scheduler: str = "slurm"


@dataclass
class BuildStrategy:
    """
    Complete build strategy for an application.
    
    Generated automatically based on application requirements
    and system capabilities.
    """
    # Selected build type
    build_type: BuildType = BuildType.CPU_RELEASE
    execution_model: ExecutionModel = ExecutionModel.CPU_ONLY
    
    # CMake configuration
    cmake_flags: Dict[str, str] = field(default_factory=dict)
    cmake_defines: List[str] = field(default_factory=list)
    
    # Compiler selection
    c_compiler: Optional[str] = None
    cxx_compiler: Optional[str] = None
    cuda_compiler: Optional[str] = None
    hip_compiler: Optional[str] = None
    fortran_compiler: Optional[str] = None
    
    # MPI configuration
    mpi_wrapper: Optional[str] = None  # mpicc, mpicxx
    
    # Required modules to load
    modules: List[str] = field(default_factory=list)
    
    # Environment variables
    env_vars: Dict[str, str] = field(default_factory=dict)
    
    # Build parallelism
    parallel_jobs: int = 4


class ApplicationAnalyzer:
    """
    Analyzes application source code to determine requirements.
    
    Scans source files for:
    - CUDA/HIP/SYCL includes and kernels
    - MPI function calls
    - OpenMP pragmas
    - Library dependencies
    """
    
    # Pattern definitions for detecting requirements
    PATTERNS = {
        'cuda': [
            r'#include\s*[<"]cuda',
            r'#include\s*[<"]cublas',
            r'#include\s*[<"]cufft',
            r'#include\s*[<"]cudnn',
            r'__global__\s+void',
            r'__device__\s+',
            r'cudaMalloc',
            r'cudaMemcpy',
            r'<<<.*>>>'
        ],
        'hip': [
            r'#include\s*[<"]hip/',
            r'#include\s*[<"]hipblas',
            r'hipMalloc',
            r'hipMemcpy',
            r'hipLaunchKernelGGL'
        ],
        'sycl': [
            r'#include\s*[<"]CL/sycl',
            r'#include\s*[<"]sycl/',
            r'sycl::queue',
            r'cl::sycl::'
        ],
        'mpi': [
            r'#include\s*[<"]mpi\.h[>"]',
            r'MPI_Init',
            r'MPI_Comm_rank',
            r'MPI_Send',
            r'MPI_Recv',
            r'MPI_Bcast'
        ],
        'openmp': [
            r'#pragma\s+omp',
            r'#include\s*[<"]omp\.h[>"]',
            r'omp_get_thread_num',
            r'omp_set_num_threads'
        ],
        'openacc': [
            r'#pragma\s+acc',
            r'acc_init',
            r'acc_get_device'
        ],
        'hdf5': [
            r'#include\s*[<"]hdf5\.h[>"]',
            r'#include\s*[<"]H5',
            r'H5Fcreate',
            r'H5Fopen'
        ],
        'fftw': [
            r'#include\s*[<"]fftw3',
            r'fftw_plan',
            r'fftw_execute'
        ],
        'petsc': [
            r'#include\s*[<"]petsc',
            r'PetscInitialize',
            r'PETSC_COMM_WORLD'
        ]
    }
    
    # File extension patterns
    CUDA_EXTENSIONS = ['.cu', '.cuh']
    HIP_EXTENSIONS = ['.hip', '.hip.cpp']
    SOURCE_EXTENSIONS = ['.c', '.cpp', '.cxx', '.cc', '.f', '.f90', '.F90', '.f95']
    
    def __init__(self, source_dir: Path):
        """
        Initialize analyzer with source directory.
        
        Args:
            source_dir: Path to application source code
        """
        self.source_dir = Path(source_dir)
        if not self.source_dir.exists():
            raise ValueError(f"Source directory does not exist: {source_dir}")
    
    def analyze(self) -> ApplicationRequirements:
        """
        Analyze source code to determine requirements.
        
        Returns:
            ApplicationRequirements with detected features
        """
        logger.info(f"Analyzing application source: {self.source_dir}")
        
        req = ApplicationRequirements()
        
        # Detect build system
        req.build_system = self._detect_build_system()
        logger.info(f"  Build system: {req.build_system}")
        
        # Collect all source files
        all_files = self._collect_source_files()
        
        # Categorize files
        for file_path in all_files:
            ext = file_path.suffix.lower()
            if ext in self.CUDA_EXTENSIONS:
                req.cuda_files.append(str(file_path))
            elif ext in self.HIP_EXTENSIONS or '.hip' in file_path.name.lower():
                req.hip_files.append(str(file_path))
            elif ext in self.SOURCE_EXTENSIONS:
                req.source_files.append(str(file_path))
        
        logger.info(f"  Source files: {len(req.source_files)}")
        logger.info(f"  CUDA files: {len(req.cuda_files)}")
        logger.info(f"  HIP files: {len(req.hip_files)}")
        
        # Scan files for patterns
        content_cache = {}
        for file_path in all_files:
            try:
                with open(file_path, 'r', errors='ignore') as f:
                    content = f.read()
                    content_cache[str(file_path)] = content
            except Exception as e:
                logger.debug(f"Could not read {file_path}: {e}")
        
        # Detect requirements from content
        all_content = '\n'.join(content_cache.values())
        
        req.requires_cuda = self._match_patterns(all_content, 'cuda') or len(req.cuda_files) > 0
        req.requires_hip = self._match_patterns(all_content, 'hip') or len(req.hip_files) > 0
        req.requires_sycl = self._match_patterns(all_content, 'sycl')
        req.requires_mpi = self._match_patterns(all_content, 'mpi')
        req.requires_openmp = self._match_patterns(all_content, 'openmp')
        req.requires_openacc = self._match_patterns(all_content, 'openacc')
        req.requires_hdf5 = self._match_patterns(all_content, 'hdf5')
        req.requires_fftw = self._match_patterns(all_content, 'fftw')
        req.requires_petsc = self._match_patterns(all_content, 'petsc')
        
        # Also check CMakeLists.txt for dependencies
        cmake_deps = self._analyze_cmake_dependencies()
        if cmake_deps.get('cuda'):
            req.requires_cuda = True
        if cmake_deps.get('mpi'):
            req.requires_mpi = True
        if cmake_deps.get('openmp'):
            req.requires_openmp = True
        
        # Determine execution model
        req.execution_model = self._determine_execution_model(req)
        
        # Calculate confidence
        req.confidence = self._calculate_confidence(req)
        
        logger.info(f"  Requirements detected:")
        logger.info(f"    MPI: {req.requires_mpi}")
        logger.info(f"    OpenMP: {req.requires_openmp}")
        logger.info(f"    CUDA: {req.requires_cuda}")
        logger.info(f"    HIP: {req.requires_hip}")
        logger.info(f"    OpenACC: {req.requires_openacc}")
        logger.info(f"  Execution model: {req.execution_model.value}")
        logger.info(f"  Confidence: {req.confidence:.2f}")
        
        return req
    
    def _detect_build_system(self) -> str:
        """Detect the build system used by the application."""
        if (self.source_dir / "CMakeLists.txt").exists():
            return "cmake"
        elif (self.source_dir / "configure.ac").exists() or (self.source_dir / "configure").exists():
            return "autotools"
        elif (self.source_dir / "Makefile").exists() or (self.source_dir / "makefile").exists():
            return "make"
        elif (self.source_dir / "meson.build").exists():
            return "meson"
        return "cmake"  # Default assumption
    
    def _collect_source_files(self) -> List[Path]:
        """Collect all source files from the source directory."""
        files = []
        extensions = self.CUDA_EXTENSIONS + self.HIP_EXTENSIONS + self.SOURCE_EXTENSIONS + ['.h', '.hpp', '.hxx']
        
        for ext in extensions:
            files.extend(self.source_dir.rglob(f"*{ext}"))
        
        # Also collect headers
        files.extend(self.source_dir.rglob("*.h"))
        files.extend(self.source_dir.rglob("*.hpp"))
        
        return files
    
    def _match_patterns(self, content: str, pattern_key: str) -> bool:
        """Check if content matches any pattern in the pattern group."""
        patterns = self.PATTERNS.get(pattern_key, [])
        for pattern in patterns:
            if re.search(pattern, content, re.IGNORECASE):
                return True
        return False
    
    def _analyze_cmake_dependencies(self) -> Dict[str, bool]:
        """Analyze CMakeLists.txt for dependencies."""
        deps = {'cuda': False, 'mpi': False, 'openmp': False, 'hip': False}
        
        cmake_file = self.source_dir / "CMakeLists.txt"
        if not cmake_file.exists():
            return deps
        
        try:
            content = cmake_file.read_text()
            
            # Check for find_package or enable_language calls
            if re.search(r'find_package\s*\(\s*CUDA', content, re.IGNORECASE):
                deps['cuda'] = True
            if re.search(r'enable_language\s*\(\s*CUDA', content, re.IGNORECASE):
                deps['cuda'] = True
            if re.search(r'find_package\s*\(\s*MPI', content, re.IGNORECASE):
                deps['mpi'] = True
            if re.search(r'find_package\s*\(\s*OpenMP', content, re.IGNORECASE):
                deps['openmp'] = True
            if re.search(r'find_package\s*\(\s*HIP', content, re.IGNORECASE):
                deps['hip'] = True
                
        except Exception as e:
            logger.debug(f"Could not analyze CMakeLists.txt: {e}")
        
        return deps
    
    def _determine_execution_model(self, req: ApplicationRequirements) -> ExecutionModel:
        """Determine the execution model based on requirements."""
        has_gpu = req.requires_cuda or req.requires_hip or req.requires_sycl or req.requires_openacc
        has_cpu_parallel = req.requires_mpi or req.requires_openmp
        
        if has_gpu and has_cpu_parallel:
            return ExecutionModel.HYBRID
        elif has_gpu:
            return ExecutionModel.GPU_ONLY
        elif has_cpu_parallel:
            return ExecutionModel.CPU_ONLY
        else:
            # Pure serial or couldn't detect - assume CPU
            return ExecutionModel.CPU_ONLY
    
    def _calculate_confidence(self, req: ApplicationRequirements) -> float:
        """Calculate confidence in the detection."""
        score = 0.5  # Base score
        
        # Increase confidence based on evidence
        if len(req.cuda_files) > 0:
            score += 0.2
        if len(req.source_files) > 0:
            score += 0.1
        if req.build_system == "cmake":
            score += 0.1
        if req.requires_mpi:
            score += 0.1
            
        return min(score, 1.0)


class SystemDetector:
    """
    Detects system hardware capabilities.
    
    Uses SLURM queries and system commands to detect:
    - CPU configuration
    - GPU configuration (NVIDIA, AMD, Intel)
    - Available compilers and modules
    """
    
    def __init__(self, partition: Optional[str] = None):
        """
        Initialize system detector.
        
        Args:
            partition: Optional SLURM partition to query
        """
        self.partition = partition
    
    def detect(self) -> SystemCapabilities:
        """
        Detect system capabilities.
        
        Returns:
            SystemCapabilities with detected hardware info
        """
        logger.info("Detecting system capabilities...")
        
        caps = SystemCapabilities()
        caps.partition = self.partition
        
        # Detect CPU configuration
        self._detect_cpu(caps)
        
        # Detect GPU configuration
        self._detect_gpu(caps)
        
        # Detect available compilers
        self._detect_compilers(caps)
        
        logger.info(f"  CPU cores/node: {caps.cpu_cores_per_node}")
        logger.info(f"  GPUs/node: {caps.gpus_per_node}")
        if caps.gpu_vendor:
            logger.info(f"  GPU vendor: {caps.gpu_vendor}")
            logger.info(f"  GPU model: {caps.gpu_model}")
        logger.info(f"  Compilers: {caps.available_compilers}")
        
        return caps
    
    def _detect_cpu(self, caps: SystemCapabilities):
        """Detect CPU configuration."""
        # Try SLURM first
        if self.partition:
            try:
                result = subprocess.run(
                    ['sinfo', '-p', self.partition, '-o', '%c', '-h'],
                    capture_output=True, text=True, timeout=10
                )
                if result.returncode == 0:
                    cores = result.stdout.strip().split('\n')[0]
                    caps.cpu_cores_per_node = int(cores)
                    return
            except Exception as e:
                logger.debug(f"SLURM CPU detection failed: {e}")
        
        # Fallback to system detection
        try:
            import multiprocessing
            caps.cpu_cores_per_node = multiprocessing.cpu_count()
        except Exception:
            caps.cpu_cores_per_node = None
    
    def _detect_gpu(self, caps: SystemCapabilities):
        """Detect GPU configuration."""
        # Try SLURM GRES first
        if self.partition:
            try:
                result = subprocess.run(
                    ['sinfo', '-p', self.partition, '-o', '%G', '-h'],
                    capture_output=True, text=True, timeout=10
                )
                if result.returncode == 0:
                    gres = result.stdout.strip().split('\n')[0]
                    # Parse gres like "gpu:4" or "gpu:a100:4"
                    match = re.search(r'gpu[:\w]*:(\d+)', gres)
                    if match:
                        caps.gpus_per_node = int(match.group(1))
            except Exception as e:
                logger.debug(f"SLURM GPU detection failed: {e}")
        
        # Detect GPU vendor and model
        self._detect_nvidia_gpu(caps)
        if caps.gpus_per_node == 0:
            self._detect_amd_gpu(caps)
        if caps.gpus_per_node == 0:
            self._detect_intel_gpu(caps)
    
    def _detect_nvidia_gpu(self, caps: SystemCapabilities):
        """Detect NVIDIA GPUs using nvidia-smi."""
        try:
            result = subprocess.run(
                ['nvidia-smi', '--query-gpu=count,name,memory.total', '--format=csv,noheader'],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                if lines and lines[0]:
                    parts = lines[0].split(',')
                    if len(parts) >= 2:
                        caps.gpus_per_node = len(lines)  # One line per GPU
                        caps.gpu_vendor = 'nvidia'
                        caps.gpu_model = parts[1].strip()
                        if len(parts) >= 3:
                            mem_str = parts[2].strip()
                            match = re.search(r'(\d+)', mem_str)
                            if match:
                                caps.gpu_memory_gb = int(match.group(1)) / 1024
        except Exception as e:
            logger.debug(f"NVIDIA GPU detection failed: {e}")
    
    def _detect_amd_gpu(self, caps: SystemCapabilities):
        """Detect AMD GPUs using rocm-smi."""
        try:
            result = subprocess.run(
                ['rocm-smi', '--showproductname'],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0:
                caps.gpu_vendor = 'amd'
                lines = result.stdout.strip().split('\n')
                # Count GPU lines
                gpu_count = sum(1 for line in lines if 'GPU' in line)
                caps.gpus_per_node = gpu_count if gpu_count > 0 else 1
                # Try to extract model
                for line in lines:
                    if 'Card series' in line:
                        caps.gpu_model = line.split(':')[-1].strip()
                        break
        except Exception as e:
            logger.debug(f"AMD GPU detection failed: {e}")
    
    def _detect_intel_gpu(self, caps: SystemCapabilities):
        """Detect Intel GPUs."""
        try:
            # Check for Intel GPU using sycl-ls or clinfo
            result = subprocess.run(
                ['sycl-ls'],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0 and 'Intel' in result.stdout:
                caps.gpu_vendor = 'intel'
                caps.gpus_per_node = result.stdout.count('Intel')
                caps.gpu_model = 'Intel GPU'
        except Exception:
            pass
    
    def _detect_compilers(self, caps: SystemCapabilities):
        """Detect available compilers."""
        compilers = {
            'gcc': ['gcc', 'g++'],
            'intel': ['icc', 'icpc', 'icx', 'icpx'],
            'clang': ['clang', 'clang++'],
            'nvcc': ['nvcc'],
            'hipcc': ['hipcc'],
            'nvc++': ['nvc++'],
        }
        
        for compiler_family, executables in compilers.items():
            for exe in executables:
                if shutil.which(exe):
                    if compiler_family not in caps.available_compilers:
                        caps.available_compilers.append(compiler_family)
                    if compiler_family in ['nvcc', 'nvc++']:
                        caps.available_cuda_compilers.append(exe)
                    if compiler_family == 'hipcc':
                        caps.available_hip_compilers.append(exe)
                    break


class BuildStrategySelector:
    """
    Selects the optimal build strategy based on application requirements
    and system capabilities.
    """
    
    def __init__(self, app_requirements: ApplicationRequirements, 
                 system_capabilities: SystemCapabilities):
        """
        Initialize strategy selector.
        
        Args:
            app_requirements: Detected application requirements
            system_capabilities: Detected system capabilities
        """
        self.app_req = app_requirements
        self.sys_caps = system_capabilities
    
    def select(self) -> BuildStrategy:
        """
        Select the optimal build strategy.
        
        Returns:
            BuildStrategy with complete build configuration
        """
        logger.info("Selecting build strategy...")
        
        strategy = BuildStrategy()
        
        # Determine build type based on requirements and capabilities
        if self.app_req.requires_cuda and self.sys_caps.gpu_vendor == 'nvidia':
            if self.app_req.requires_mpi:
                strategy.build_type = BuildType.HYBRID_CUDA
                strategy.execution_model = ExecutionModel.HYBRID
            else:
                strategy.build_type = BuildType.GPU_CUDA
                strategy.execution_model = ExecutionModel.GPU_ONLY
        elif self.app_req.requires_hip and self.sys_caps.gpu_vendor == 'amd':
            if self.app_req.requires_mpi:
                strategy.build_type = BuildType.HYBRID_HIP
                strategy.execution_model = ExecutionModel.HYBRID
            else:
                strategy.build_type = BuildType.GPU_HIP
                strategy.execution_model = ExecutionModel.GPU_ONLY
        else:
            strategy.build_type = BuildType.CPU_RELEASE
            strategy.execution_model = ExecutionModel.CPU_ONLY
        
        # Configure compilers
        self._configure_compilers(strategy)
        
        # Configure CMake flags
        self._configure_cmake_flags(strategy)
        
        # Configure modules
        self._configure_modules(strategy)
        
        # Configure environment
        self._configure_environment(strategy)
        
        logger.info(f"  Build type: {strategy.build_type.value}")
        logger.info(f"  Execution model: {strategy.execution_model.value}")
        logger.info(f"  Compilers: C={strategy.c_compiler}, CXX={strategy.cxx_compiler}")
        if strategy.cuda_compiler:
            logger.info(f"  CUDA compiler: {strategy.cuda_compiler}")
        
        return strategy
    
    def _configure_compilers(self, strategy: BuildStrategy):
        """Configure compiler selection."""
        # C/C++ compilers
        if 'intel' in self.sys_caps.available_compilers:
            strategy.c_compiler = 'icc'
            strategy.cxx_compiler = 'icpc'
        elif 'gcc' in self.sys_caps.available_compilers:
            strategy.c_compiler = 'gcc'
            strategy.cxx_compiler = 'g++'
        elif 'clang' in self.sys_caps.available_compilers:
            strategy.c_compiler = 'clang'
            strategy.cxx_compiler = 'clang++'
        
        # CUDA compiler
        if strategy.build_type in [BuildType.GPU_CUDA, BuildType.HYBRID_CUDA]:
            if 'nvc++' in self.sys_caps.available_cuda_compilers:
                strategy.cuda_compiler = 'nvc++'
            elif 'nvcc' in self.sys_caps.available_cuda_compilers:
                strategy.cuda_compiler = 'nvcc'
        
        # HIP compiler
        if strategy.build_type in [BuildType.GPU_HIP, BuildType.HYBRID_HIP]:
            if self.sys_caps.available_hip_compilers:
                strategy.hip_compiler = self.sys_caps.available_hip_compilers[0]
        
        # MPI wrappers
        if self.app_req.requires_mpi:
            strategy.mpi_wrapper = 'mpicc'
    
    def _configure_cmake_flags(self, strategy: BuildStrategy):
        """Configure CMake flags based on build type."""
        flags = {
            'CMAKE_BUILD_TYPE': 'Release'
        }
        
        # Compiler settings
        if strategy.c_compiler:
            flags['CMAKE_C_COMPILER'] = strategy.c_compiler
        if strategy.cxx_compiler:
            flags['CMAKE_CXX_COMPILER'] = strategy.cxx_compiler
        
        # CUDA settings
        if strategy.build_type in [BuildType.GPU_CUDA, BuildType.HYBRID_CUDA]:
            flags['CMAKE_CUDA_COMPILER'] = strategy.cuda_compiler or 'nvcc'
            
            # Detect CUDA architecture
            arch = self._detect_cuda_architecture()
            if arch:
                flags['CMAKE_CUDA_ARCHITECTURES'] = arch
        
        # HIP settings
        if strategy.build_type in [BuildType.GPU_HIP, BuildType.HYBRID_HIP]:
            flags['CMAKE_HIP_COMPILER'] = strategy.hip_compiler or 'hipcc'
        
        # OpenMP
        if self.app_req.requires_openmp:
            flags['CMAKE_CXX_FLAGS'] = '-fopenmp'
        
        strategy.cmake_flags = flags
    
    def _detect_cuda_architecture(self) -> Optional[str]:
        """Detect CUDA architecture from GPU."""
        try:
            result = subprocess.run(
                ['nvidia-smi', '--query-gpu=compute_cap', '--format=csv,noheader'],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode == 0:
                cap = result.stdout.strip().split('\n')[0]
                # Convert "8.0" to "80"
                arch = cap.replace('.', '')
                return arch
        except Exception:
            pass
        return None
    
    def _configure_modules(self, strategy: BuildStrategy):
        """Configure required modules."""
        modules = []
        
        # Base compiler modules
        if 'intel' in self.sys_caps.available_compilers:
            modules.append('intel')
        elif 'gcc' in self.sys_caps.available_compilers:
            modules.append('gcc')
        
        # CUDA modules
        if strategy.build_type in [BuildType.GPU_CUDA, BuildType.HYBRID_CUDA]:
            if 'nvc++' in self.sys_caps.available_cuda_compilers:
                modules.append('nvhpc')
            else:
                modules.append('cuda')
        
        # HIP modules
        if strategy.build_type in [BuildType.GPU_HIP, BuildType.HYBRID_HIP]:
            modules.append('rocm')
        
        # MPI
        if self.app_req.requires_mpi:
            modules.append('openmpi')
        
        # Libraries
        if self.app_req.requires_hdf5:
            modules.append('hdf5')
        if self.app_req.requires_fftw:
            modules.append('fftw')
        
        strategy.modules = modules
    
    def _configure_environment(self, strategy: BuildStrategy):
        """Configure environment variables."""
        env = {}
        
        # OpenMP threads
        if self.app_req.requires_openmp:
            env['OMP_NUM_THREADS'] = '1'  # Will be set properly at runtime
        
        # CUDA
        if strategy.build_type in [BuildType.GPU_CUDA, BuildType.HYBRID_CUDA]:
            if 'CUDA_HOME' in os.environ:
                env['CUDA_HOME'] = os.environ['CUDA_HOME']
        
        strategy.env_vars = env


def auto_configure_build(source_dir: Path, 
                         partition: Optional[str] = None,
                         hardware_type: Optional[str] = None) -> Tuple[BuildStrategy, ApplicationRequirements, SystemCapabilities]:
    """
    Automatically configure build strategy for an application.
    
    This is the main entry point for automated build configuration.
    
    Args:
        source_dir: Path to application source code
        partition: Optional SLURM partition for hardware detection
        hardware_type: Optional hint ('cpu' or 'gpu') to override detection
        
    Returns:
        Tuple of (BuildStrategy, ApplicationRequirements, SystemCapabilities)
    """
    logger.info("=" * 60)
    logger.info("AUTOMATED BUILD CONFIGURATION")
    logger.info("=" * 60)
    
    # Analyze application
    analyzer = ApplicationAnalyzer(source_dir)
    app_req = analyzer.analyze()
    
    # Detect system capabilities
    detector = SystemDetector(partition)
    sys_caps = detector.detect()
    
    # Override execution model if hardware_type is specified
    if hardware_type:
        if hardware_type.lower() == 'gpu':
            if sys_caps.gpus_per_node > 0:
                app_req.execution_model = ExecutionModel.GPU_ONLY
            else:
                logger.warning(f"GPU requested but no GPUs detected on partition {partition}")
        elif hardware_type.lower() == 'cpu':
            app_req.execution_model = ExecutionModel.CPU_ONLY
    
    # Select build strategy
    selector = BuildStrategySelector(app_req, sys_caps)
    strategy = selector.select()
    
    logger.info("=" * 60)
    logger.info("BUILD CONFIGURATION COMPLETE")
    logger.info("=" * 60)
    
    return strategy, app_req, sys_caps


# Convenience function for integration with existing code
def get_build_flags_for_hardware(source_dir: Path,
                                  partition: Optional[str] = None,
                                  hardware_type: str = 'cpu') -> Dict[str, str]:
    """
    Get CMake flags for the detected hardware configuration.
    
    Args:
        source_dir: Path to source code
        partition: SLURM partition
        hardware_type: 'cpu' or 'gpu'
        
    Returns:
        Dictionary of CMake flags
    """
    strategy, _, _ = auto_configure_build(source_dir, partition, hardware_type)
    return strategy.cmake_flags


if __name__ == '__main__':
    import sys
    logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
    
    if len(sys.argv) < 2:
        print("Usage: python build_strategy.py <source_dir> [partition] [hardware_type]")
        sys.exit(1)
    
    source_dir = Path(sys.argv[1])
    partition = sys.argv[2] if len(sys.argv) > 2 else None
    hardware_type = sys.argv[3] if len(sys.argv) > 3 else None
    
    strategy, app_req, sys_caps = auto_configure_build(source_dir, partition, hardware_type)
    
    print("\n" + "=" * 60)
    print("RESULT SUMMARY")
    print("=" * 60)
    print(f"Application: {source_dir.name}")
    print(f"Execution Model: {strategy.execution_model.value}")
    print(f"Build Type: {strategy.build_type.value}")
    print(f"\nCMake Flags:")
    for key, value in strategy.cmake_flags.items():
        print(f"  -D{key}={value}")
    print(f"\nModules: {', '.join(strategy.modules)}")
