"""
Build System Discovery - REFINED VERSION

Analyzes source code to automatically detect:
- Build system type (CMake, Make, Autotools, etc.)
- Required modules and compilers
- Build flags and commands

CRITICAL FIXES:
- Removed hardcoded compiler/module patterns.
- Added support for a wider range of HPC compilers (Intel, AMD, PGI/NVHPC).
- Modularized detection logic for better maintainability.
"""

import re
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

# --- Configuration for Flexible Detection ---
# This configuration would ideally be loaded from a system-wide YAML/JSON file
# to allow system administrators to easily update patterns without modifying code.
# For this fix, we define it internally.

COMPILER_PATTERNS = {
    'cuda': r'(?:cuda|CUDA)[\s/_-]?(\d+[._]\d+(?:[._]\d+)?)',
    'gcc': r'(?:gcc|GCC)[\s/_-]?(\d+[._]\d+(?:[._]\d+)?)',
    'cmake': r'(?:cmake|CMake)[\s/_-]?(\d+[._]\d+(?:[._]\d+)?)',
    'intel': r'intel|icc|icpc|ifort',
    'amd': r'amd|aocc|hip|rocm',
    'nvhpc': r'nvhpc|nvc|nvc\+\+|nvfortran|pgi',
    'openmpi': r'openmpi',
    'mpich': r'mpich',
}

MODULE_PATTERNS = [
    r'module load ([a-zA-Z0-9_\-\.\/]+)',
    r'module\s+load\s+([a-zA-Z0-9_\-\.\/]+)',
]

BUILD_SYSTEM_FILES = {
    "cmake": ["CMakeLists.txt"],
    "autotools": ["configure.ac", "configure.in"],
    "make": ["Makefile", "makefile"],
    "spack": ["spack.yaml", "package.py"],
    "easybuild": ["*.eb"],
}

# --- Data Structures ---
@dataclass
class BuildInfo:
    """Information about build system and requirements."""
    build_system: str  # 'cmake', 'make', 'autotools', 'spack', 'easybuild', 'custom'
    required_modules: List[str] = field(default_factory=list)
    compiler_requirements: Dict[str, str] = field(default_factory=dict)  # e.g., {'cuda': '11.8', 'gcc': '11.3'}
    mpi_requirement: Optional[str] = None  # e.g., 'openmpi', 'mpich'
    build_flags: Dict[str, str] = field(default_factory=dict)
    build_commands: List[str] = field(default_factory=list)
    source_files: List[str] = field(default_factory=list)
    readme_content: str = ""
    detected_features: List[str] = field(default_factory=list)


# --- BuildDiscovery Class ---
class BuildDiscovery:
    """Discovers build system requirements from source code."""
    
    def __init__(self, system_config: Optional[Dict] = None):
        """
        Initialize BuildDiscovery.
        
        Args:
            system_config: Optional system configuration dict with GPU info
        """
        self.system_config = system_config or {}
        self._gpu_arch_cache = None
    
    def analyze_source(self, source_dir: Path) -> BuildInfo:
        """
        Analyze source directory to detect build requirements.
        
        Args:
            source_dir: Path to source code directory
            
        Returns:
            BuildInfo object with detected information
        """
        source_dir = Path(source_dir)
        
        logger.info(f"Analyzing source directory: {source_dir}")
        
        # 1. Read README
        readme_content = self._read_readme(source_dir)
        
        # 2. Detect build system
        build_system = self._detect_build_system(source_dir)
        logger.info(f"  Build system: {build_system}")
        
        # 3. Extract information from README and source files
        required_modules = self._extract_modules(readme_content)
        compiler_requirements = self._extract_compiler_requirements(readme_content)
        mpi_requirement = self._extract_mpi_requirement(readme_content)
        build_flags = self._extract_build_flags(readme_content, build_system, source_dir)
        build_commands = self._extract_build_commands(readme_content, build_system)
        source_files = self._scan_source_files(source_dir)
        detected_features = self._detect_features(readme_content, source_files)
        
        build_info = BuildInfo(
            build_system=build_system,
            required_modules=required_modules,
            compiler_requirements=compiler_requirements,
            mpi_requirement=mpi_requirement,
            build_flags=build_flags,
            build_commands=build_commands,
            readme_content=readme_content,
            source_files=source_files,
            detected_features=detected_features
        )
        
        logger.info(f"  Required modules: {len(required_modules)}")
        logger.info(f"  Build flags: {len(build_flags)}")
        
        return build_info
    
    # --- Modularized Detection Methods ---
    
    def _detect_build_system(self, source_dir: Path) -> str:
        """Detect the build system used by the project based on file presence."""
        for system, files in BUILD_SYSTEM_FILES.items():
            for file_pattern in files:
                if '*' in file_pattern:
                    if list(source_dir.glob(file_pattern)):
                        return system
                elif (source_dir / file_pattern).exists():
                    return system
        return "custom"

    def _read_readme(self, source_dir: Path) -> str:
        """Read README file if it exists."""
        readme_patterns = ["README.md", "README.rst", "README.txt", "README", "readme.md"]
        
        for pattern in readme_patterns:
            readme_path = source_dir / pattern
            if readme_path.exists():
                try:
                    return readme_path.read_text(encoding='utf-8', errors='ignore')
                except Exception as e:
                    logger.warning(f"Could not read {readme_path}: {e}")
        return ""

    def _scan_source_files(self, source_dir: Path) -> List[str]:
        """Scan for source files to detect languages and features."""
        extensions = {'.c', '.cpp', '.cc', '.cxx', '.cu', '.f', '.f90', '.f95', '.py', '.hip', '.cl', '.sycl'}
        source_files = []
        
        try:
            for ext in extensions:
                source_files.extend([str(f.relative_to(source_dir)) 
                                   for f in source_dir.rglob(f"*{ext}")])
        except Exception as e:
            logger.warning(f"Error scanning source files: {e}")
        
        return source_files[:100]

    def _extract_modules(self, readme_content: str) -> List[str]:
        """Extract required modules from README using configurable patterns."""
        modules = []
        
        for pattern in MODULE_PATTERNS:
            matches = re.findall(pattern, readme_content, re.MULTILINE)
            modules.extend(matches)
        
        # Remove duplicates while preserving order
        seen = set()
        unique_modules = []
        for mod in modules:
            if mod not in seen:
                seen.add(mod)
                unique_modules.append(mod)
        
        return unique_modules

    def _extract_compiler_requirements(self, readme_content: str) -> Dict[str, str]:
        """Extract compiler version requirements from README using configurable patterns."""
        requirements = {}
        
        for compiler, pattern in COMPILER_PATTERNS.items():
            match = re.search(pattern, readme_content, re.IGNORECASE)
            if match:
                # For patterns with a capture group (e.g., version), extract the group
                if match.groups():
                    requirements[compiler] = match.group(1)
                # For patterns without a capture group (e.g., just a name), set a generic requirement
                else:
                    requirements[compiler] = "required"
        
        return requirements

    def _extract_mpi_requirement(self, readme_content: str) -> Optional[str]:
        """Extract MPI implementation requirement."""
        # Check for specific MPI implementations (prioritize specific over generic)
        if re.search(r'\bcray.?mpich\b', readme_content, re.IGNORECASE):
            return 'cray-mpich'
        elif re.search(r'\bintel.?mpi\b', readme_content, re.IGNORECASE):
            return 'intel-mpi'
        elif re.search(r'\bopenmpi\b', readme_content, re.IGNORECASE):
            return 'openmpi'
        elif re.search(r'\bmpich\b', readme_content, re.IGNORECASE):
            return 'mpich'
        elif re.search(r'\bmpi\b', readme_content, re.IGNORECASE):
            return 'mpi'  # Generic MPI
        
        return None

    def _detect_features(self, readme_content: str, source_files: List[str]) -> List[str]:
        """Detect parallel programming features (CUDA, OpenMP, OpenACC, HIP, SYCL)."""
        features = []
        
        # Check source files
        if any(f.endswith('.cu') for f in source_files):
            features.append('cuda')
        if any(f.endswith('.hip') for f in source_files):
            features.append('hip')
        if any(f.endswith('.sycl') for f in source_files):
            features.append('sycl')
        
        # Check README content
        content = readme_content.lower()
        if 'openmp' in content:
            features.append('openmp')
        if 'openacc' in content:
            features.append('openacc')
        if 'cuda' in content and 'cuda' not in features:
            features.append('cuda')
        if 'hip' in content and 'hip' not in features:
            features.append('hip')
        if 'rocm' in content and 'hip' not in features:
            features.append('hip') # Assume ROCm implies HIP for build purposes
        if 'sycl' in content and 'sycl' not in features:
            features.append('sycl')
            
        return sorted(list(set(features)))

    def _get_gpu_architecture_from_system(self) -> Optional[str]:
        """
        Get GPU architecture from system configuration.
        
        Returns:
            GPU architecture string (e.g., "80" for sm_80) or None
        """
        if self._gpu_arch_cache is not None:
            return self._gpu_arch_cache
        
        try:
            if self.system_config:
                systems = self.system_config.get('systems', [])
                for system in systems:
                    partitions = system.get('partitions', [])
                    for partition in partitions:
                        devices = partition.get('devices', [])
                        for device in devices:
                            if device.get('type') == 'gpu':
                                arch = device.get('arch', '')
                                # Extract numeric part from "sm_80" -> "80"
                                if arch.startswith('sm_'):
                                    self._gpu_arch_cache = arch[3:]
                                    logger.info(f"Detected GPU architecture from system config: sm_{self._gpu_arch_cache}")
                                    return self._gpu_arch_cache
        except Exception as e:
            logger.debug(f"Could not extract GPU arch from system config: {e}")
        
        self._gpu_arch_cache = None
        return None
    
    def _detect_cuda_architecture_from_source(self, source_dir: Path) -> Optional[str]:
        """
        Detect CUDA architecture from CMakeLists.txt or source files.
        
        Returns:
            Architecture string (e.g., "70;80") or None
        """
        cmake_file = source_dir / "CMakeLists.txt"
        if cmake_file.exists():
            try:
                content = cmake_file.read_text(encoding='utf-8', errors='ignore')
                
                # Look for CMAKE_CUDA_ARCHITECTURES
                arch_pattern = r'set\s*\(\s*CMAKE_CUDA_ARCHITECTURES\s+["\']?([0-9;]+)["\']?\s*\)'
                match = re.search(arch_pattern, content, re.IGNORECASE)
                if match:
                    arch = match.group(1)
                    logger.info(f"Found CUDA architecture in CMakeLists.txt: {arch}")
                    return arch
            except Exception as e:
                logger.debug(f"Error reading CMakeLists.txt: {e}")
        
        return None

    def _extract_build_flags(self, readme_content: str, build_system: str, source_dir: Path) -> Dict[str, str]:
        """
        Extract build flags from README and inject dynamic flags (e.g., CUDA arch).
        
        Args:
            readme_content: README file content
            build_system: Detected build system
            source_dir: Source directory path
            
        Returns:
            Dictionary of build flags
        """
        flags = {}
        
        if build_system == 'cmake':
            # Extract CMake flags from README
            cmake_pattern = r'-D([A-Z_]+)=([^\s]+)'
            matches = re.findall(cmake_pattern, readme_content)
            for key, value in matches:
                flags[key] = value
        
        elif build_system == 'make':
            # Extract Make variables from README
            make_pattern = r'make\s+([A-Z_]+)=([^\s]+)'
            matches = re.findall(make_pattern, readme_content)
            for key, value in matches:
                flags[key] = value
        
        # --- Dynamic Flag Injection (CUDA Architecture) ---
        
        # Check if CUDA is required based on features or compiler requirements
        cuda_required = 'cuda' in self._detect_features(readme_content, self._scan_source_files(source_dir))
        
        if cuda_required:
            flags['ENABLE_CUDA'] = 'ON'
            
            # CRITICAL FIX: Do NOT hardcode CUDA architectures! (Copied from previous fix)
            cuda_arch = None
            
            # 1. Try to detect from source
            source_arch = self._detect_cuda_architecture_from_source(source_dir)
            if source_arch:
                cuda_arch = source_arch
                logger.info(f"Using CUDA architecture from source: {cuda_arch}")
            
            # 2. Try to get from system configuration
            if not cuda_arch:
                system_arch = self._get_gpu_architecture_from_system()
                if system_arch:
                    cuda_arch = system_arch
                    logger.info(f"Using CUDA architecture from system config: {cuda_arch}")
            
            # 3. Use "native" for CMake 3.24+ (auto-detection)
            # This check must happen *after* checking source and system config
            if not cuda_arch:
                cmake_version = self._extract_compiler_requirements(readme_content).get('cmake', '')
                if cmake_version:
                    try:
                        major, minor = map(int, cmake_version.split('.')[:2])
                        if major > 3 or (major == 3 and minor >= 24):
                            cuda_arch = "native"
                            logger.info("Using CMAKE_CUDA_ARCHITECTURES=native (auto-detection)")
                    except:
                        pass
            
            # 4. Last resort: use "all-major" for maximum compatibility
            if not cuda_arch:
                cuda_arch = "all-major"
                logger.warning("No CUDA architecture detected, using 'all-major' for maximum compatibility")
            
            flags['CMAKE_CUDA_ARCHITECTURES'] = cuda_arch
        
        return flags
    
    def _extract_build_commands(self, readme_content: str, build_system: str) -> List[str]:
        """
        Extract build commands from README.
        
        Args:
            readme_content: README file content
            build_system: Detected build system
            
        Returns:
            List of build commands
        """
        commands = []
        
         # Look for code blocks with build commands
        code_block_pattern = r'```[a-zA-Z]*\s*(.*?)\s*```'
        code_blocks = re.findall(code_block_pattern, readme_content, re.DOTALL)
        
        for block in code_blocks:
            lines = block.strip().split('\n')
            for line in lines:
                line = line.strip()
                # Filter for build-related commands, ignoring comments and module loads
                if line and not line.startswith('#') and not line.startswith('module load'):
                    if any(cmd in line for cmd in ['cmake', 'make', 'configure', 'build', 'spack', 'eb']):
                        commands.append(line)
        
        return commands


def get_build_recommendations(build_info: BuildInfo) -> Dict[str, Any]:
    """
    Generate build recommendations from BuildInfo.
    
    Args:
        build_info: BuildInfo object
        
    Returns:
        Dictionary with build recommendations
    """
    recommendations = {
        'build_system': build_info.build_system,
        'modules': build_info.required_modules,
        'cmake_flags': build_info.build_flags,
        'parallel_jobs': 8,
        'env_vars': {}
    }
    
    # Add MPI compiler wrappers if MPI is required
    if build_info.mpi_requirement:
        recommendations['env_vars']['CC'] = 'mpicc'
        recommendations['env_vars']['CXX'] = 'mpicxx'
        recommendations['env_vars']['FC'] = 'mpifort'
    
    # Add CUDA paths if CUDA is required
    if 'cuda' in build_info.detected_features:
        recommendations['env_vars']['CUDA_HOME'] = '${CUDA_HOME}'
        if build_info.build_system == 'cmake' and 'CMAKE_CUDA_COMPILER' not in build_info.build_flags:
            recommendations['cmake_flags']['CMAKE_CUDA_COMPILER'] = 'nvcc'
    
    # Add OpenMP flags if required
    if 'openmp' in build_info.detected_features:
        # Only add if not already present in build flags
        if 'CMAKE_C_FLAGS' not in recommendations['cmake_flags']:
            recommendations['cmake_flags']['CMAKE_C_FLAGS'] = '-fopenmp'
        if 'CMAKE_CXX_FLAGS' not in recommendations['cmake_flags']:
            recommendations['cmake_flags']['CMAKE_CXX_FLAGS'] = '-fopenmp'
    
    return recommendations
