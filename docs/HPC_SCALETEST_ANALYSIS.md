# HPC-ScaleTest Framework Analysis and Improvement Recommendations

## Executive Summary

This document presents a comprehensive review of the HPC-ScaleTest framework, identifying limitations, assumptions, and architectural improvements to achieve the stated goal of being **system, compiler, and application agnostic**.

---

## 1. Architectural Issues

### 1.1 Plugin Registry Design Flaws

**Current State:**
```python
class PluginRegistry:
    _directives: Dict[str, Callable] = {}  # Class-level mutable state
    _features: Dict[str, Type] = {}
    _backends: Dict[str, Type] = {}
    ...
```

**Problems:**
1. Class-level mutable dictionaries persist across tests, causing state pollution
2. No mechanism to unregister plugins (testing cleanup)
3. No plugin versioning or dependency management
4. Build systems not included in registry despite having an interface

**Recommended Fix:**
```python
class PluginRegistry:
    """Singleton registry with instance-level state and lifecycle management."""
    
    _instance: Optional['PluginRegistry'] = None
    
    def __init__(self):
        self._plugins: Dict[str, Dict[str, PluginEntry]] = {
            'scheduler': {},
            'launcher': {},
            'build_system': {},
            'module_system': {},
            'vcs': {},  # NEW: Version control systems
            'result_parser': {},  # NEW: Custom parsers
        }
        self._dependencies: Dict[str, List[str]] = {}
    
    @classmethod
    def instance(cls) -> 'PluginRegistry':
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance
    
    @classmethod
    def reset(cls):
        """Reset for testing isolation."""
        cls._instance = None
    
    def register(self, category: str, name: str, cls: Type, 
                 version: str = "1.0.0", depends: List[str] = None):
        """Register with version and dependency tracking."""
        ...
    
    def unregister(self, category: str, name: str):
        """Remove plugin for testing."""
        ...
```

### 1.2 Factory Pattern Inconsistency

**Current State:** Mix of registry lookups and hardcoded imports:
```python
# Some use registry
if has_launcher(backend):
    return get_launcher(backend, options)

# Others use hardcoded imports
if backend == BuildBackend.CMAKE:
    from backends.builds.cmake import CMakeBackend
    return CMakeBackend(options)
```

**Problems:**
1. Cannot add custom schedulers without modifying factory code
2. Build systems bypass registry entirely
3. Inconsistent extension mechanism

**Recommended Fix:** Unified registry-based factory:
```python
class BackendFactory:
    @staticmethod
    def create(category: str, backend: str, options: Dict = None) -> Any:
        """Create any backend type from unified registry."""
        registry = PluginRegistry.instance()
        plugin_class = registry.get(category, backend)
        if plugin_class is None:
            raise ValueError(f"Unknown {category}: {backend}")
        return plugin_class(options)
```

### 1.3 Abstract Interface Gaps

**Missing Methods in `BuildSystemInterface`:**
```python
class BuildSystemInterface(ABC):
    # EXISTING: configure, build, install, clean
    
    # MISSING - Should be added:
    @abstractmethod
    def find_executable(self, build_dir: Path, source_dir: Path = None) -> Optional[Path]:
        """Auto-detect compiled binary."""
        pass
    
    @abstractmethod
    def build_and_find(self, source_dir: Path, build_dir: Path, 
                       flags: Dict = None) -> Optional[Path]:
        """Complete build workflow with binary detection."""
        pass
    
    @abstractmethod
    def get_build_info(self, build_dir: Path) -> Dict:
        """Return build metadata (compiler, flags, dependencies)."""
        pass
```

**Missing Interface: `VersionControlInterface`:**
```python
class VersionControlInterface(ABC):
    """Abstract interface for version control systems."""
    
    @abstractmethod
    def clone(self, url: str, target: Path, **options) -> Path:
        """Clone repository."""
        pass
    
    @abstractmethod
    def checkout(self, target: Path, ref: str) -> bool:
        """Checkout specific branch/tag/commit."""
        pass
    
    @abstractmethod
    def get_revision(self, target: Path) -> str:
        """Get current revision identifier."""
        pass
    
    @abstractmethod
    def update(self, target: Path) -> bool:
        """Update to latest."""
        pass
```

---

## 2. Generality Issues

### 2.1 Scaling Engine Tied to PIC Simulations

**Current State:** Parameter names are hardcoded for Particle-in-Cell codes:
```python
self.param_names = [
    'Lx', 'Ly', 'Lz',           # iPIC3D specific
    'nxc', 'nyc', 'nzc',        # iPIC3D specific
    'XLEN', 'YLEN', 'ZLEN',     # iPIC3D specific
    'npcelx', 'npcely', 'npcelz' # iPIC3D specific
]
```

**Impact:** Cannot scale OpenFOAM, LAMMPS, GROMACS, or any other application without code changes.

**Recommended Fix:** Fully parameterized scaling engine:
```yaml
# User-defined variable mapping in YAML
variable_mapping:
  domain_size:
    x: "Lx"         # or "xdim" for another code
    y: "Ly"
    z: "Lz"
  grid_cells:
    x: "nxc"        # or "nx" or "IMAX"
    y: "nyc"
    z: "nzc"
  mpi_decomposition:
    x: "XLEN"       # or "npx" or "NPROCS_X"
    y: "YLEN"
    z: "ZLEN"
  
  # Custom parameters for specific codes
  custom:
    timesteps: "nsteps"
    output_freq: "iout"
```

```python
class GenericScalingEngine:
    """Application-agnostic scaling engine."""
    
    def __init__(self, variable_mapping: Dict):
        self.var_map = variable_mapping
    
    def scale(self, scaling_type: str, nodes: List[int], config: Dict) -> List[Dict]:
        """Generate scaled configurations using user-defined variable names."""
        scaled_configs = []
        
        for node_count in nodes:
            if scaling_type == "weak":
                scaled = self._weak_scale(node_count, config)
            elif scaling_type == "strong":
                scaled = self._strong_scale(node_count, config)
            scaled_configs.append(scaled)
        
        return scaled_configs
```

### 2.2 README Analyzer Limitations

**Current State:** Simple regex patterns that miss many cases:
```python
MPI_PATTERNS = [
    r'\b(mpi|openmpi|mpich|intelmpi)\b',
    ...
]
```

**Problems:**
1. Cannot parse CI/CD config files (.gitlab-ci.yml, .github/workflows)
2. Cannot detect Spack/EasyBuild recipes
3. Cannot analyze Dockerfiles/Singularity definitions
4. No NLP/ML for understanding free-form text

**Recommended Fix:** Multi-source analyzer:
```python
class ProjectAnalyzer:
    """Analyze project from multiple sources."""
    
    def __init__(self, source_dir: Path):
        self.analyzers = [
            ReadmeAnalyzer(source_dir),
            CMakeAnalyzer(source_dir),
            CIConfigAnalyzer(source_dir),      # NEW
            ContainerFileAnalyzer(source_dir), # NEW
            SpackRecipeAnalyzer(source_dir),   # NEW
            EasyBuildRecipeAnalyzer(source_dir), # NEW
        ]
    
    def analyze(self) -> BuildInfo:
        """Combine results from all analyzers with confidence weighting."""
        results = []
        for analyzer in self.analyzers:
            try:
                result = analyzer.analyze()
                results.append(result)
            except Exception as e:
                logger.debug(f"{analyzer.__class__.__name__} failed: {e}")
        
        return self._merge_results(results)
```

### 2.3 Input File Format Assumptions

**Current State:** Input parser assumes specific delimiters:
```python
patterns = [
    rf'{re.escape(param_name)}\s*[=:]\s*([^\s#;]+)',  # = or : delimiter
    rf'{re.escape(param_name)}\s+(\S+)',              # space-separated
]
```

**Problems:** Cannot handle:
- XML configuration files
- JSON/YAML configuration
- Fortran namelists (`&namelist ... /`)
- HDF5/NetCDF attributes
- Binary configuration files

**Recommended Fix:** Format-specific parsers:
```python
class InputParserFactory:
    """Factory for format-specific parsers."""
    
    parsers = {
        '.inp': KeyValueParser,
        '.yaml': YAMLParser,
        '.yml': YAMLParser,
        '.json': JSONParser,
        '.xml': XMLParser,
        '.nml': NamelistParser,     # Fortran namelist
        '.toml': TOMLParser,
    }
    
    @classmethod
    def get_parser(cls, file_path: Path) -> InputParser:
        suffix = file_path.suffix.lower()
        parser_class = cls.parsers.get(suffix, KeyValueParser)
        return parser_class()
```

---

## 3. Missing Capabilities

### 3.1 Version Control System Support

**Current State:** Only basic Git clone:
```python
def _clone_repository(self, git_url: str) -> Path:
    subprocess.run(["git", "clone", git_url, str(target_dir)])
```

**Missing:**
- Branch/tag/commit checkout
- Shallow cloning
- Submodule support
- Authentication (SSH keys, tokens)
- Mercurial, SVN support
- Sparse checkout

**Recommended Implementation:**
```python
class GitBackend(VersionControlInterface):
    def clone(self, url: str, target: Path, **options) -> Path:
        cmd = ["git", "clone"]
        
        if options.get('shallow', False):
            cmd.extend(["--depth", "1"])
        
        if options.get('branch'):
            cmd.extend(["--branch", options['branch']])
        
        if options.get('recursive', True):
            cmd.append("--recursive")
        
        cmd.extend([url, str(target)])
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise VCSError(f"Clone failed: {result.stderr}")
        
        # Handle specific commit checkout
        if options.get('commit'):
            self.checkout(target, options['commit'])
        
        return target
    
    def checkout(self, target: Path, ref: str) -> bool:
        """Checkout branch, tag, or commit."""
        result = subprocess.run(
            ["git", "checkout", ref],
            cwd=target,
            capture_output=True
        )
        return result.returncode == 0
```

### 3.2 Container Support

**Missing:** No Singularity/Docker integration for reproducible builds.

**Recommended Implementation:**
```python
class ContainerBuildBackend(BuildSystemInterface):
    """Build inside containers for reproducibility."""
    
    def __init__(self, options: Dict = None):
        self.runtime = options.get('runtime', 'singularity')  # or 'docker'
        self.image = options.get('image')
    
    def build_in_container(self, source_dir: Path, build_dir: Path, 
                           flags: Dict = None) -> Optional[Path]:
        """Execute build inside container."""
        if self.runtime == 'singularity':
            cmd = [
                "singularity", "exec",
                "--bind", f"{source_dir}:/src",
                "--bind", f"{build_dir}:/build",
                self.image,
                "/bin/bash", "-c",
                "cd /build && cmake /src && make -j"
            ]
        else:
            cmd = [
                "docker", "run", "--rm",
                "-v", f"{source_dir}:/src",
                "-v", f"{build_dir}:/build",
                "-w", "/build",
                self.image,
                "bash", "-c", "cmake /src && make -j"
            ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        return build_dir if result.returncode == 0 else None
```

### 3.3 Additional Scheduler Support

**Current State:** Only SLURM and LOCAL implemented.

**Missing Schedulers:**
- PBS/Torque
- SGE/UGE
- LSF
- Flux

**Recommended:** Abstract scheduler directives:
```python
class PBSScheduler(SchedulerInterface):
    """PBS/Torque scheduler backend."""
    
    DIRECTIVE_MAP = {
        'nodes': '#PBS -l nodes=',
        'time': '#PBS -l walltime=',
        'queue': '#PBS -q ',
        'account': '#PBS -A ',
        'job_name': '#PBS -N ',
        'output': '#PBS -o ',
        'error': '#PBS -e ',
    }
    
    def generate_job_script(self, job_config, resource_config, 
                            command, env_setup) -> str:
        script = "#!/bin/bash\n"
        script += f"{self.DIRECTIVE_MAP['nodes']}{job_config.num_nodes}:ppn={resource_config.procs_per_node}\n"
        script += f"{self.DIRECTIVE_MAP['time']}{resource_config.time_limit}\n"
        # ... etc
        return script
```

### 3.4 Result Parser Extensibility

**Current State:** Hardcoded regex patterns for timing extraction:
```python
time_patterns = [
    r'(?:Time|Wall\s+time|Runtime|Elapsed):\s*(\d+(?:\.\d+)?)',
    r'Total time: (\d+(?:\.\d+)?) seconds',
]
```

**Recommended:** User-defined extraction patterns:
```yaml
# In run.yaml
result_parsing:
  timing:
    pattern: "Total simulation time: ([\d.]+) seconds"
    group: 1
    unit: "seconds"
  
  performance:
    pattern: "Performance: ([\d.]+) GFLOPS"
    group: 1
    unit: "GFLOPS"
  
  memory:
    pattern: "Peak memory: ([\d.]+) GB"
    group: 1
    unit: "GB"
  
  custom_metrics:
    - name: "particles_per_second"
      pattern: "Throughput: ([\d.e+]+) particles/s"
      group: 1
```

```python
class ConfigurableResultParser(ResultParserInterface):
    def __init__(self, patterns: Dict):
        self.patterns = patterns
    
    def parse_output(self, output_file: Path) -> Dict:
        content = output_file.read_text()
        results = {}
        
        for metric_name, config in self.patterns.items():
            match = re.search(config['pattern'], content)
            if match:
                results[metric_name] = {
                    'value': float(match.group(config.get('group', 1))),
                    'unit': config.get('unit', '')
                }
        
        return results
```

### 3.5 Performance Regression Detection

**Missing:** No capability to compare results across runs.

**Recommended Implementation:**
```python
class PerformanceAnalyzer:
    """Detect performance regressions."""
    
    def __init__(self, baseline_results: Path):
        self.baseline = self._load_baseline(baseline_results)
    
    def analyze(self, current_results: Dict, threshold: float = 0.1) -> Dict:
        """Compare current results to baseline."""
        analysis = {
            'regressions': [],
            'improvements': [],
            'unchanged': []
        }
        
        for metric, current_value in current_results.items():
            if metric not in self.baseline:
                continue
            
            baseline_value = self.baseline[metric]
            change = (current_value - baseline_value) / baseline_value
            
            if change > threshold:
                analysis['regressions'].append({
                    'metric': metric,
                    'baseline': baseline_value,
                    'current': current_value,
                    'change_pct': change * 100
                })
            elif change < -threshold:
                analysis['improvements'].append({...})
            else:
                analysis['unchanged'].append({...})
        
        return analysis
```

---

## 4. Robustness Issues

### 4.1 Error Handling and Recovery

**Current State:** Many operations fail silently or with generic errors.

**Recommended:** Structured error handling with recovery:
```python
class HPCScaleTestError(Exception):
    """Base exception for all framework errors."""
    pass

class BuildError(HPCScaleTestError):
    """Build-related errors."""
    def __init__(self, message: str, build_log: str = None, 
                 suggestions: List[str] = None):
        super().__init__(message)
        self.build_log = build_log
        self.suggestions = suggestions or []

class ModuleError(HPCScaleTestError):
    """Module system errors."""
    def __init__(self, message: str, missing_modules: List[str] = None,
                 available_alternatives: List[str] = None):
        super().__init__(message)
        self.missing_modules = missing_modules
        self.available_alternatives = available_alternatives
```

### 4.2 Module Validation

**Current State:** Modules are loaded without verification.

**Recommended:** Pre-flight module validation:
```python
class ModuleValidator:
    """Validate modules before build/run."""
    
    def validate_modules(self, modules: List[str]) -> Tuple[bool, List[str], Dict]:
        """Check if all modules are available."""
        available = []
        missing = []
        alternatives = {}
        
        for module in modules:
            if self.module_system.is_module_available(module):
                available.append(module)
            else:
                missing.append(module)
                # Find alternatives
                base_name = module.split('/')[0]
                alts = self.module_system.list_available_modules(f"{base_name}/*")
                if alts:
                    alternatives[module] = alts
        
        return len(missing) == 0, missing, alternatives
```

### 4.3 Checkpoint/Restart Support

**Missing:** No support for long-running tests that may be interrupted.

**Recommended Implementation:**
```python
class CheckpointManager:
    """Manage test checkpoints for restart capability."""
    
    def __init__(self, run_dir: Path):
        self.checkpoint_file = run_dir / ".checkpoint.json"
    
    def save_checkpoint(self, state: Dict):
        """Save current state."""
        with open(self.checkpoint_file, 'w') as f:
            json.dump({
                'timestamp': datetime.now().isoformat(),
                'state': state
            }, f)
    
    def load_checkpoint(self) -> Optional[Dict]:
        """Load previous state if exists."""
        if self.checkpoint_file.exists():
            with open(self.checkpoint_file) as f:
                return json.load(f)['state']
        return None
    
    def can_resume(self) -> bool:
        """Check if a valid checkpoint exists."""
        return self.checkpoint_file.exists()
```

---

## 5. Recommended Architectural Changes

### 5.1 Proposed Package Structure

```
hpc_scaletest/
├── core/
│   ├── __init__.py
│   ├── interfaces/           # All abstract interfaces
│   │   ├── scheduler.py
│   │   ├── launcher.py
│   │   ├── build_system.py
│   │   ├── vcs.py            # NEW
│   │   ├── result_parser.py
│   │   └── container.py      # NEW
│   ├── registry.py           # Unified plugin registry
│   ├── factory.py            # Registry-based factory
│   ├── config.py
│   └── exceptions.py         # Structured exceptions
│
├── plugins/                   # All plugins in one location
│   ├── schedulers/
│   │   ├── slurm.py
│   │   ├── pbs.py
│   │   ├── sge.py
│   │   └── local.py
│   ├── launchers/
│   │   ├── srun.py
│   │   ├── mpirun.py
│   │   └── jsrun.py
│   ├── build_systems/
│   │   ├── cmake.py
│   │   ├── make.py
│   │   ├── autotools.py
│   │   ├── meson.py          # NEW
│   │   └── bazel.py          # NEW
│   ├── vcs/                   # NEW
│   │   ├── git.py
│   │   ├── mercurial.py
│   │   └── svn.py
│   ├── modules/
│   │   ├── lmod.py
│   │   ├── tmod.py
│   │   └── spack.py          # NEW
│   └── containers/           # NEW
│       ├── singularity.py
│       └── docker.py
│
├── engine/
│   ├── orchestrator.py
│   ├── runner.py
│   ├── scaling/              # Refactored scaling
│   │   ├── base.py           # Generic scaling engine
│   │   ├── weak.py
│   │   └── strong.py
│   └── analysis/
│       ├── result_parser.py
│       └── regression.py     # NEW
│
├── utils/
│   ├── analyzers/            # Refactored analyzers
│   │   ├── readme.py
│   │   ├── cmake.py
│   │   ├── ci_config.py      # NEW
│   │   └── container.py      # NEW
│   ├── parsers/              # Format-specific parsers
│   │   ├── keyvalue.py
│   │   ├── namelist.py
│   │   ├── yaml_parser.py
│   │   └── xml_parser.py
│   └── validators/
│
└── cli/                       # Command-line interface
    ├── __init__.py
    └── main.py
```

### 5.2 Plugin Auto-Discovery

```python
# core/registry.py
import importlib
import pkgutil

def autodiscover_plugins(package_name: str = 'hpc_scaletest.plugins'):
    """Automatically discover and register all plugins."""
    package = importlib.import_module(package_name)
    
    for importer, modname, ispkg in pkgutil.walk_packages(
        path=package.__path__,
        prefix=package.__name__ + '.',
        onerror=lambda x: None
    ):
        try:
            importlib.import_module(modname)
            logger.debug(f"Loaded plugin module: {modname}")
        except Exception as e:
            logger.warning(f"Failed to load plugin {modname}: {e}")
```

### 5.3 Configuration Schema Validation

```python
from pydantic import BaseModel, validator
from typing import Optional, List, Dict

class ScalingConfig(BaseModel):
    type: str  # "weak" or "strong"
    max_nodes: int
    initial_cells: Optional[List[int]] = None
    initial_domain: Optional[List[float]] = None
    
    @validator('type')
    def validate_type(cls, v):
        if v not in ['weak', 'strong']:
            raise ValueError('scaling type must be "weak" or "strong"')
        return v

class HPCScaleTestConfig(BaseModel):
    """Validated configuration schema."""
    repository: str
    build_system: str = "cmake"
    scaling: ScalingConfig
    modules: List[str] = []
    
    # Variable mapping for application-agnostic scaling
    variable_mapping: Optional[Dict[str, Dict[str, str]]] = None
```

---

## 6. Summary of Recommendations

### High Priority (Core Functionality)
1. **Refactor Plugin Registry** - Instance-based with lifecycle management
2. **Add BuildSystemInterface.find_executable** - Align interface with implementation
3. **Create VersionControlInterface** - Abstract VCS operations
4. **Implement GenericScalingEngine** - Remove hardcoded parameter names

### Medium Priority (Extensibility)
5. **Add PBS/SGE/LSF Schedulers** - Broader HPC support
6. **Implement Container Build Backend** - Reproducibility
7. **Create Format-Specific Input Parsers** - Handle diverse config formats
8. **Add Custom Result Parser Support** - User-defined metrics

### Lower Priority (Robustness)
9. **Structured Exception Hierarchy** - Better error handling
10. **Checkpoint/Restart Support** - Long-running test recovery
11. **Performance Regression Detection** - Automated comparison
12. **Multi-source Project Analyzer** - CI configs, containers, recipes

### Documentation
13. **Plugin Development Guide** - How to add custom plugins
14. **Variable Mapping Reference** - Application-specific examples
15. **Troubleshooting Guide** - Common issues and solutions

---

## 7. Conclusion

The HPC-ScaleTest framework provides a solid foundation for automated scaling studies but requires architectural improvements to achieve true system, compiler, and application agnosticism. The key changes are:

1. **Plugin Registry Redesign** - Enable true extensibility
2. **Variable Mapping System** - Remove application-specific hardcoding
3. **Multi-format Support** - Handle diverse input file formats
4. **Broader Scheduler Support** - Beyond SLURM
5. **Container Integration** - Reproducible builds

These changes would transform HPC-ScaleTest from a PIC-simulation-focused tool into a truly generic HPC benchmarking framework.
