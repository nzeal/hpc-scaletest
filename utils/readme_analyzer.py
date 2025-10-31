"""
README analyzer for automatic dependency and build system detection.
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Set
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class BuildInfo:
    """Information extracted from README and project files."""
    build_system: Optional[str] = None
    modules: List[str] = field(default_factory=list)
    dependencies: List[str] = field(default_factory=list)
    compilers: List[str] = field(default_factory=list)
    mpi_required: bool = False
    cuda_required: bool = False
    openmp_required: bool = False
    build_commands: List[str] = field(default_factory=list)
    cmake_flags: Dict[str, str] = field(default_factory=dict)
    confidence_score: float = 0.0


class ReadmeAnalyzer:
    """Analyzes README and project files to extract build information."""
    
    # Common patterns for dependency detection
    COMPILER_PATTERNS = {
        'gcc': r'\b(gcc|g\+\+|gfortran)\b',
        'intel': r'\b(icc|icpc|ifort|intel)\b',
        'clang': r'\b(clang|clang\+\+)\b',
        'nvcc': r'\b(nvcc|cuda)\b'
    }
    
    MPI_PATTERNS = [
        r'\b(mpi|openmpi|mpich|intelmpi|mpicc|mpicxx|mpirun|mpiexec)\b',
        r'\bMPI(_|\s)',
        r'#include\s*<mpi\.h>'
    ]
    
    CUDA_PATTERNS = [
        r'\b(cuda|nvcc|cudart|cublas|cufft)\b',
        r'\bCUDA(_|\s)',
        r'#include\s*<cuda'
    ]
    
    OPENMP_PATTERNS = [
        r'\b(openmp|omp)\b',
        r'\bOpenMP',
        r'#pragma\s+omp'
    ]
    
    MODULE_PATTERNS = {
        'gcc': [r'gcc/[\d.]+', r'gnu/[\d.]+'],
        'intel': [r'intel/[\d.]+', r'intel-oneapi/[\d.]+'],
        'openmpi': [r'openmpi/[\d.]+'],
        'mpich': [r'mpich/[\d.]+'],
        'cuda': [r'cuda/[\d.]+', r'cudatoolkit/[\d.]+'],
        'hdf5': [r'hdf5/[\d.]+'],
        'fftw': [r'fftw/[\d.]+'],
        'boost': [r'boost/[\d.]+']
    }
    
    def __init__(self, source_dir: Path):
        """
        Initialize README analyzer.
        
        Args:
            source_dir: Path to source code directory
        """
        self.source_dir = source_dir
        self.readme_files = self._find_readme_files()
        self.build_files = self._find_build_files()
    
    def analyze(self) -> BuildInfo:
        """
        Analyze README and build files to extract build information.
        
        Returns:
            BuildInfo object with detected information
        """
        info = BuildInfo()
        
        # Detect build system
        info.build_system = self._detect_build_system()
        
        # Analyze README files
        readme_content = self._read_all_readmes()
        if readme_content:
            info.mpi_required = self._detect_mpi(readme_content)
            info.cuda_required = self._detect_cuda(readme_content)
            info.openmp_required = self._detect_openmp(readme_content)
            info.compilers = self._detect_compilers(readme_content)
            info.modules = self._detect_modules(readme_content)
            info.build_commands = self._extract_build_commands(readme_content)
            info.dependencies = self._extract_dependencies(readme_content)
        
        # Analyze build files for additional information
        if info.build_system == "cmake":
            info.cmake_flags = self._analyze_cmake_files()
        
        # Calculate confidence score
        info.confidence_score = self._calculate_confidence(info)
        
        logger.info(f"Build system detected: {info.build_system} (confidence: {info.confidence_score:.2f})")
        logger.info(f"MPI required: {info.mpi_required}, CUDA required: {info.cuda_required}")
        logger.info(f"Detected modules: {info.modules}")
        
        return info
    
    def _find_readme_files(self) -> List[Path]:
        """Find README files in source directory."""
        readme_patterns = ["README*", "readme*", "READ.ME", "INSTALL*"]
        readme_files = []
        for pattern in readme_patterns:
            readme_files.extend(self.source_dir.glob(pattern))
        return readme_files
    
    def _find_build_files(self) -> Dict[str, List[Path]]:
        """Find build system files."""
        build_files = {
            'cmake': list(self.source_dir.glob("**/CMakeLists.txt")),
            'makefile': list(self.source_dir.glob("**/[Mm]akefile*")),
            'autotools': list(self.source_dir.glob("**/configure*")),
            'meson': list(self.source_dir.glob("**/meson.build"))
        }
        return build_files
    
    def _read_all_readmes(self) -> str:
        """Read and concatenate all README files."""
        content = []
        for readme in self.readme_files:
            try:
                content.append(readme.read_text(errors='ignore'))
            except Exception as e:
                logger.warning(f"Could not read {readme}: {e}")
        return "\n".join(content)
    
    def _detect_build_system(self) -> Optional[str]:
        """Detect build system from project files."""
        # Priority order
        if self.build_files['cmake']:
            return "cmake"
        elif self.build_files['autotools']:
            return "autotools"
        elif self.build_files['makefile']:
            return "make"
        elif self.build_files['meson']:
            return "meson"
        return None
    
    def _detect_mpi(self, content: str) -> bool:
        """Detect if MPI is required."""
        content_lower = content.lower()
        for pattern in self.MPI_PATTERNS:
            if re.search(pattern, content_lower, re.IGNORECASE):
                return True
        return False
    
    def _detect_cuda(self, content: str) -> bool:
        """Detect if CUDA is required."""
        content_lower = content.lower()
        for pattern in self.CUDA_PATTERNS:
            if re.search(pattern, content_lower, re.IGNORECASE):
                return True
        return False
    
    def _detect_openmp(self, content: str) -> bool:
        """Detect if OpenMP is required."""
        content_lower = content.lower()
        for pattern in self.OPENMP_PATTERNS:
            if re.search(pattern, content_lower, re.IGNORECASE):
                return True
        return False
    
    def _detect_compilers(self, content: str) -> List[str]:
        """Detect required compilers."""
        compilers = []
        content_lower = content.lower()
        for compiler, pattern in self.COMPILER_PATTERNS.items():
            if re.search(pattern, content_lower, re.IGNORECASE):
                compilers.append(compiler)
        return compilers
    
    def _detect_modules(self, content: str) -> List[str]:
        """Detect module names from README."""
        modules = []
        content_lower = content.lower()
        
        for module_type, patterns in self.MODULE_PATTERNS.items():
            for pattern in patterns:
                matches = re.findall(pattern, content_lower, re.IGNORECASE)
                if matches:
                    # Use the first match found
                    modules.append(matches[0])
                    break
        
        # Also look for explicit module load commands
        module_load_pattern = r'module\s+load\s+([\w/.-]+)'
        explicit_modules = re.findall(module_load_pattern, content_lower, re.IGNORECASE)
        modules.extend(explicit_modules)
        
        # Remove duplicates and filter
        unique_modules = list(set(modules))
        
        # Filter out modules that are likely to not exist (very specific versions)
        # Keep only base module names or common versions
        filtered_modules = []
        for module in unique_modules:
            # Skip modules with very specific patch versions (e.g., x.y.z where z > 10)
            parts = module.split('/')
            if len(parts) == 2:
                name, version = parts
                version_parts = version.split('.')
                try:
                    # Keep module if version seems reasonable (not too specific)
                    if len(version_parts) <= 2 or (len(version_parts) == 3 and int(version_parts[2]) <= 10):
                        filtered_modules.append(module)
                    else:
                        logger.debug(f"Filtering out overly specific module: {module}")
                except (ValueError, IndexError):
                    # If version can't be parsed, keep the module
                    filtered_modules.append(module)
            else:
                filtered_modules.append(module)
        
        return filtered_modules
    
    def _extract_build_commands(self, content: str) -> List[str]:
        """Extract build commands from README."""
        commands = []
        
        # Look for code blocks with build commands
        cmake_patterns = [
            r'cmake\s+[^\n]+',
            r'mkdir\s+build\s*&&\s*cd\s+build',
            r'make\s+[^\n]*',
            r'\./configure\s+[^\n]+'
        ]
        
        for pattern in cmake_patterns:
            matches = re.findall(pattern, content, re.IGNORECASE)
            commands.extend(matches)
        
        return commands
    
    def _extract_dependencies(self, content: str) -> List[str]:
        """Extract dependency names from README."""
        dependencies = []
        
        # Common dependency indicators
        dep_patterns = [
            r'(?:requires?|depends?\s+on|needs?)\s*:?\s*([^\n]+)',
            r'(?:dependencies|prerequisites)\s*:?\s*([^\n]+)'
        ]
        
        for pattern in dep_patterns:
            matches = re.findall(pattern, content, re.IGNORECASE)
            for match in matches:
                # Split on common delimiters
                deps = re.split(r'[,;]', match)
                dependencies.extend([d.strip() for d in deps if d.strip()])
        
        return dependencies
    
    def _analyze_cmake_files(self) -> Dict[str, str]:
        """Analyze CMakeLists.txt for flags and options."""
        flags = {}
        
        for cmake_file in self.build_files['cmake']:
            try:
                content = cmake_file.read_text(errors='ignore')
                
                # Look for find_package calls
                packages = re.findall(r'find_package\s*\(\s*(\w+)', content, re.IGNORECASE)
                for pkg in packages:
                    if pkg.upper() == 'MPI':
                        flags['CMAKE_C_COMPILER'] = 'mpicc'
                        flags['CMAKE_CXX_COMPILER'] = 'mpicxx'
                    elif pkg.upper() == 'CUDA':
                        flags['CUDA_TOOLKIT_ROOT_DIR'] = '${CUDA_HOME}'
                
                # Look for options
                options = re.findall(r'option\s*\(\s*(\w+)\s+([^\)]+)\)', content, re.IGNORECASE)
                for opt_name, opt_desc in options:
                    if 'USE' in opt_name.upper() or 'ENABLE' in opt_name.upper():
                        logger.debug(f"Found CMake option: {opt_name} - {opt_desc}")
                
            except Exception as e:
                logger.warning(f"Could not analyze {cmake_file}: {e}")
        
        return flags
    
    def _calculate_confidence(self, info: BuildInfo) -> float:
        """Calculate confidence score for detected information."""
        score = 0.0
        
        # Build system detection (30%)
        if info.build_system:
            score += 0.3
        
        # Module detection (20%)
        if info.modules:
            score += 0.2
        
        # Compiler detection (15%)
        if info.compilers:
            score += 0.15
        
        # Build commands (15%)
        if info.build_commands:
            score += 0.15
        
        # Dependency detection (20%)
        if info.mpi_required or info.cuda_required or info.openmp_required:
            score += 0.2
        
        return min(score, 1.0)
    
    def generate_build_recommendation(self, info: BuildInfo) -> Dict[str, any]:
        """
        Generate build recommendations based on analyzed information.
        
        Args:
            info: BuildInfo object
            
        Returns:
            Dictionary with build recommendations
        """
        recommendations = {
            'build_system': info.build_system or 'cmake',
            'modules': info.modules or [],
            'cmake_flags': info.cmake_flags or {},
            'parallel_jobs': 8,
            'env_vars': {}
        }
        
        # Add compiler environment variables if MPI is required
        if info.mpi_required:
            recommendations['env_vars']['CC'] = 'mpicc'
            recommendations['env_vars']['CXX'] = 'mpicxx'
            recommendations['env_vars']['FC'] = 'mpifort'
        
        # Add CUDA paths if CUDA is required
        if info.cuda_required:
            recommendations['env_vars']['CUDA_HOME'] = '${CUDA_HOME}'
            recommendations['cmake_flags']['CMAKE_CUDA_COMPILER'] = 'nvcc'
        
        # Add OpenMP flags if required
        if info.openmp_required:
            recommendations['cmake_flags']['CMAKE_C_FLAGS'] = '-fopenmp'
            recommendations['cmake_flags']['CMAKE_CXX_FLAGS'] = '-fopenmp'
        
        return recommendations
