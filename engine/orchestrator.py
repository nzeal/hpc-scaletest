"""
High-level orchestrator for end-to-end HPC scaling workflow.
Automates code acquisition, analysis, compilation, testing, and reporting.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, List
from dataclasses import dataclass

from core.test_definition import Test
from core.types import ScalingType, SchedulerBackend, LauncherBackend, ModuleBackend, BuildBackend
from engine.runner import TestRunner
from utils.code_acquisition import CodeAcquisition
from utils.readme_analyzer import ReadmeAnalyzer, BuildInfo
from utils.report_generator import ReportGenerator
from utils.parameter_suggestion import ParameterSuggestion, suggest_parameters
from utils.system_info import get_partition_info

logger = logging.getLogger(__name__)


@dataclass
class OrchestratorConfig:
    """Configuration for the orchestrator."""
    # Source configuration
    source: str  # Local path or Git URL
    
    # Test configuration
    scaling_type: str = "strong"  # "strong" or "weak"
    hardware_type: str = "cpu"  # "cpu" or "gpu"
    
    # Scaling parameters
    max_nodes: int = 4
    initial_procs: tuple = (2, 2, 2)
    initial_domain: Optional[tuple] = None
    initial_cells: Optional[tuple] = None
    particles_per_cell: Optional[tuple] = None  # (npcelx, npcely, npcelz)
    
    # Scaling factor (if defined, enables full weak scaling mode)
    scaling_factor: Optional[float] = None  # e.g., 2 (doubles per step)
    
    # Scaling dimensions: 1 for 1D (X only), 2 for 2D (X→Y→X→Y), 3 for 3D (X→Y→Z→X→Y→Z)
    scaling_dimensions: int = 2  # Default to 2D scaling
    
    # Scaling factors (if defined, overrides node-based scaling)
    weak_scaling_factors: Optional[List[float]] = None  # e.g., [1, 2, 4, 8]
    strong_scaling_factors: Optional[List[float]] = None  # e.g., [1, 2, 4, 8]
    
    # GENERIC: Variable mapping for application-agnostic input parsing
    variable_map: Optional[Dict] = None  # Maps parameter types to variable names
    
    # Resource configuration
    procs_per_node: int = 128
    gpus_per_node: int = 0
    memory_per_node_gb: float = 256.0  # Memory in GB
    time_limit: str = "02:00:00"
    partition: str = "X_usr_prod"
    qos: Optional[str] = None  # QoS specification for job submission
    qos_mapping: Optional[Dict[str, Dict]] = None  # QoS threshold mapping by node count
    # Example: {"small": {"max_nodes": 16, "qos": "normal"}, "large": {"min_nodes": 17, "qos": "dcgp_qos_bprod"}}
    account: str = "cin_X"
    
    # Backend configuration (auto-detected if None)
    scheduler: str = "slurm"
    launcher: Optional[str] = None  # Auto-select based on scheduler
    module_system: str = "lmod"
    
    # Build configuration
    build_system: Optional[str] = None  # Auto-detect from README
    modules: Optional[List[str]] = None  # Auto-detect from README
    build_flags: Optional[Dict[str, str]] = None
    
    # Input file configuration
    input_file: Optional[str] = None  # Path to input file or directory
    input_file_name: Optional[str] = None  # Name of input file (e.g., "input.txt")
    auto_detect_input: bool = True  # Auto-detect input files from source
    
    # Behavior configuration
    auto_submit_jobs: bool = True
    cleanup_after_build: bool = False
    generate_reports: bool = True
    
    # Output configuration
    output_dir: Path = Path("output")
    workspace_dir: Path = Path("workspace")
    test_name: Optional[str] = None
    
    # LLM configuration (for intelligent parameter mapping)
    use_llm_discovery: bool = False
    openai_api_key: Optional[str] = None
    llm_model: str = "gpt-4"
    parameter_mapping: Optional[Dict[str, List[str]]] = None  # Custom parameter name mappings


class HPCOrchestrator:
    """
    High-level orchestrator for automated HPC scaling workflow.
    
    This class provides a simple interface for running end-to-end HPC
    scaling tests with minimal user input.
    """
    
    def __init__(self, config: OrchestratorConfig):
        """
        Initialize the orchestrator.
        
        Args:
            config: Orchestrator configuration
        """
        self.config = config
        self.source_dir: Optional[Path] = None
        self.is_cloned: bool = False
        self.build_info: Optional[BuildInfo] = None
        self.test: Optional[Test] = None
        
        # Create output and workspace directories
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        self.config.workspace_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("Initialized HPC Orchestrator")
        logger.info(f"  Source: {self.config.source}")
        logger.info(f"  Scaling: {self.config.scaling_type}")
        logger.info(f"  Hardware: {self.config.hardware_type}")
        logger.info(f"  Max nodes: {self.config.max_nodes}")
    
    def run(self) -> bool:
        """
        Execute the complete workflow.
        
        Returns:
            True if workflow completed successfully
        """
        try:
            # Step 1: Acquire source code
            logger.info("\n" + "="*60)
            logger.info("STEP 1: CODE ACQUISITION")
            logger.info("="*60)
            self._acquire_code()
            
            # Step 2: Analyze code and dependencies
            logger.info("\n" + "="*60)
            logger.info("STEP 2: CODE ANALYSIS & DEPENDENCY DETECTION")
            logger.info("="*60)
            self._analyze_code()
            
            # Step 3: Create test configuration
            logger.info("\n" + "="*60)
            logger.info("STEP 3: TEST CONFIGURATION")
            logger.info("="*60)
            self._create_test()
            
            # Step 4: Run scaling tests
            logger.info("\n" + "="*60)
            logger.info("STEP 4: SCALING TEST EXECUTION")
            logger.info("="*60)
            success = self._run_tests()
            
            if not success:
                logger.error("Test execution failed")
                return False
            
            # Step 5: Generate reports (only if jobs completed)
            if self.config.generate_reports:
                if hasattr(self, 'jobs_completed') and not self.jobs_completed:
                    logger.info("\n" + "="*60)
                    logger.info("SKIPPING STEP 5: REPORT GENERATION")
                    logger.info("="*60)
                    logger.warning("⚠ Jobs are still running or were interrupted")
                    logger.info("Reports will be generated after jobs complete")
                    logger.info("To generate reports later, run:")
                    logger.info(f"  python utils/generate_report.py {self.run_dir if hasattr(self, 'run_dir') else 'output/...'}")
                else:
                    logger.info("\n" + "="*60)
                    logger.info("STEP 5: REPORT GENERATION")
                    logger.info("="*60)
                    self._generate_reports()
            
            logger.info("\n" + "="*60)
            logger.info("WORKFLOW COMPLETED SUCCESSFULLY")
            logger.info("="*60)
            return True
            
        except Exception as e:
            logger.error(f"Workflow failed: {e}", exc_info=True)
            return False
        finally:
            # Cleanup if needed
            if self.config.cleanup_after_build and self.is_cloned:
                self._cleanup()
    
    def _acquire_code(self):
        """Step 1: Acquire source code from local path or Git repository."""
        acquisition = CodeAcquisition(self.config.workspace_dir)
        self.source_dir, self.is_cloned = acquisition.acquire(self.config.source)
        
        logger.info(f"Source code location: {self.source_dir}")
        logger.info(f"Source type: {'Git clone' if self.is_cloned else 'Local path'}")
    
    def _analyze_code(self):
        """Step 2: Analyze README and build files for dependencies."""
        analyzer = ReadmeAnalyzer(self.source_dir)
        self.build_info = analyzer.analyze()
        
        logger.info(f"Build system: {self.build_info.build_system}")
        logger.info(f"MPI required: {self.build_info.mpi_required}")
        logger.info(f"CUDA required: {self.build_info.cuda_required}")
        logger.info(f"OpenMP required: {self.build_info.openmp_required}")
        
        if self.build_info.modules:
            logger.info(f"Detected modules: {', '.join(self.build_info.modules)}")
        
        if self.build_info.compilers:
            logger.info(f"Detected compilers: {', '.join(self.build_info.compilers)}")
        
        # Generate build recommendations
        recommendations = analyzer.generate_build_recommendation(self.build_info)
        logger.info(f"Build recommendations generated (confidence: {self.build_info.confidence_score:.2f})")
        
        # Override config with detected values if not explicitly set
        if self.config.build_system is None:
            self.config.build_system = recommendations['build_system']
        
        if self.config.modules is None:
            self.config.modules = recommendations['modules']
        
        if self.config.build_flags is None:
            self.config.build_flags = recommendations['cmake_flags']
    
    def _create_test(self):
        """Step 3: Create test configuration based on analysis."""
        # Determine test name
        if self.config.test_name is None:
            self.config.test_name = self.source_dir.name
        
        # Find executable name (try common patterns)
        exe_candidates = [
            self.source_dir / "build" / self.source_dir.name,
            self.source_dir / "build" / "bin" / self.source_dir.name,
            self.source_dir / self.source_dir.name
        ]
        
        exe_path = None
        for candidate in exe_candidates:
            if candidate.exists():
                exe_path = candidate
                break
        
        # Create test instance with absolute path
        # Use resolve() to convert relative paths to absolute paths
        # so job scripts can run after changing directories
        if exe_path:
            command = [str(exe_path.resolve())]
        else:
            command = ["./a.out"]
        
        self.test = Test(
            name=self.config.test_name,
            command=command,
            output_dir=self.config.output_dir
        )
        
        # Configure backend
        launcher = self.config.launcher
        if launcher is None:
            # Auto-select launcher based on scheduler
            if self.config.scheduler == "slurm":
                launcher = "srun"
            else:
                launcher = "mpirun"
        
        self.test.set_backend(
            scheduler=self.config.scheduler,
            launcher=launcher,
            module_system=self.config.module_system
        )
        
        # Configure resources
        gpus = self.config.gpus_per_node
        if self.config.hardware_type.lower() == "gpu" and gpus == 0:
            gpus = 1  # Default to 1 GPU per node for GPU tests
        
        self.test.set_resources(
            max_nodes=self.config.max_nodes,
            procs_per_node=self.config.procs_per_node,
            gpus_per_node=gpus,
            time_limit=self.config.time_limit,
            partition=self.config.partition,
            qos=self.config.qos,
            qos_mapping=self.config.qos_mapping,
            account=self.config.account
        )
        
        # Configure scaling
        self.test.set_scaling(
            scaling_type=self.config.scaling_type,
            max_nodes=self.config.max_nodes,
            initial_procs=self.config.initial_procs,
            initial_domain=self.config.initial_domain,
            initial_cells=self.config.initial_cells
        )
        
        # Set variable mapping (for generic input parsing)
        if self.config.variable_map:
            self.test.scaling_config.variable_map = self.config.variable_map
            logger.info(f"  Variable mapping configured for generic input parsing")
        
        # Set scaling factor if defined (enables full weak scaling mode)
        if self.config.scaling_factor:
            self.test.scaling_config.scaling_factor = self.config.scaling_factor
            logger.info(f"  Scaling factor: {self.config.scaling_factor} (full weak scaling enabled)")
        
        # Set scaling dimensions (1D, 2D, or 3D scaling pattern)
        if hasattr(self.config, 'scaling_dimensions'):
            self.test.scaling_config.scaling_dimensions = self.config.scaling_dimensions
            mode_str = {1: "1D (X only)", 2: "2D (X→Y→X→Y, Z constant)", 3: "3D (X→Y→Z cycling)"}
            logger.info(f"  Scaling dimensions: {mode_str.get(self.config.scaling_dimensions, self.config.scaling_dimensions)}")
        
        # Set scaling factors if defined
        if self.config.weak_scaling_factors:
            self.test.scaling_config.weak_scaling_factors = self.config.weak_scaling_factors
            logger.info(f"  Weak scaling factors: {self.config.weak_scaling_factors}")
        
        if self.config.strong_scaling_factors:
            self.test.scaling_config.strong_scaling_factors = self.config.strong_scaling_factors
            logger.info(f"  Strong scaling factors: {self.config.strong_scaling_factors}")
        
        # Set particles per cell if specified
        if self.config.particles_per_cell:
            # Store in scaling config for use during input generation
            self.test.scaling_config.particles_per_cell = self.config.particles_per_cell
            logger.info(f"  Particles per cell: {self.config.particles_per_cell}")
        
        # Set modules
        if self.config.modules:
            self.test.set_modules(self.config.modules)
        
        # Set environment variables
        env_vars = {}
        if self.build_info and self.build_info.openmp_required:
            env_vars['OMP_NUM_THREADS'] = '1'  # Disable OpenMP for pure MPI
        
        if env_vars:
            self.test.set_env(env_vars)
        
        # Configure build
        if self.source_dir:
            self.test.build_config.source_dir = self.source_dir
            self.test.build_config.build_dir = self.source_dir / "build"
            self.test.build_config.install_dir = self.config.output_dir / "install"
            
            if self.config.build_system:
                build_system_map = {
                    'cmake': BuildBackend.CMAKE,
                    'make': BuildBackend.MAKE,
                    'autotools': BuildBackend.AUTOTOOLS
                }
                self.test.backend_config.build_system = build_system_map.get(
                    self.config.build_system.lower(),
                    BuildBackend.CMAKE
                )
            
            if self.config.build_flags:
                self.test.build_config.build_flags = self.config.build_flags
            
            # Set executable name
            self.test.build_config.executable_name = self.source_dir.name
        
        # Configure input files
        input_file_path = self._setup_input_files()
        if input_file_path:
            self.test.input_file = input_file_path
            logger.info(f"  Input file: {input_file_path}")
        
        # Set auto-submit behavior
        self.test.set_auto_submit(self.config.auto_submit_jobs)
        
        # Set LLM configuration for intelligent parameter mapping
        self.test.use_llm_discovery = self.config.use_llm_discovery
        self.test.openai_api_key = self.config.openai_api_key
        self.test.llm_model = self.config.llm_model
        self.test.parameter_mapping = self.config.parameter_mapping
        
        if self.config.use_llm_discovery:
            logger.info(f"  LLM parameter discovery: ENABLED")
            logger.info(f"  LLM model: {self.config.llm_model}")
        if self.config.parameter_mapping:
            logger.info(f"  Custom parameter mappings: {len(self.config.parameter_mapping)} types")
        
        logger.info(f"Test configured: {self.test.name}")
        logger.info(f"  Backend: {self.config.scheduler}/{launcher}")
        logger.info(f"  Resources: {self.config.max_nodes} nodes × {self.config.procs_per_node} procs/node")
        logger.info(f"  Scaling: {self.config.scaling_type}")
        logger.info(f"  Auto-submit: {self.config.auto_submit_jobs}")
        
        # Show parameter suggestions
        self._show_parameter_suggestions()
        
        # Show partition hardware info
        self._show_partition_info()
    
    def _setup_input_files(self) -> Optional[Path]:
        """Detect and prepare input files for the test."""
        # If user specified an input file, use it
        if self.config.input_file:
            input_path = Path(self.config.input_file)
            if input_path.exists():
                logger.info(f"Using user-specified input file: {input_path}")
                return input_path.resolve()
            else:
                logger.warning(f"Specified input file not found: {input_path}")
        
        # Auto-detect input files from source directory
        if self.config.auto_detect_input and self.source_dir:
            logger.info("Auto-detecting input files...")
            
            # Common input directory names
            input_dirs = ["inputfiles", "input", "inputs", "testcases", "examples", "data"]
            
            for dir_name in input_dirs:
                input_dir = self.source_dir / dir_name
                if input_dir.exists() and input_dir.is_dir():
                    logger.info(f"Found input directory: {input_dir}")
                    
                    # If user specified an input file name, look for it
                    if self.config.input_file_name:
                        input_file = input_dir / self.config.input_file_name
                        if input_file.exists():
                            logger.info(f"Found specified input file: {input_file}")
                            return input_file.resolve()
                    
                    # Otherwise, return the directory
                    logger.info(f"Using input directory: {input_dir}")
                    return input_dir.resolve()
            
            # Check for common input file names in source root
            common_input_files = [
                "input.txt", "input.dat", "input.in", "input.xml",
                "config.txt", "config.dat", "config.in", "config.xml",
                "parameters.txt", "parameters.dat", "parameters.in",
                "stdin", "os-stdin"
            ]
            
            for filename in common_input_files:
                input_file = self.source_dir / filename
                if input_file.exists():
                    logger.info(f"Found input file: {input_file}")
                    return input_file.resolve()
            
            logger.info("No input files detected (application may not need them)")
        
        return None
    
    def _show_parameter_suggestions(self):
        """Display dynamically computed parameter suggestions."""
        try:
            logger.info("\n" + "="*60)
            logger.info("PARAMETER SUGGESTIONS (based on detected hardware)")
            logger.info("="*60)
            
            # Create parameter suggester
            is_gpu_run = self.config.hardware_type.lower() == 'gpu'
            suggester = ParameterSuggestion(
                cores_per_node=self.config.procs_per_node,
                memory_gb=float(self.config.memory_per_node_gb) if hasattr(self.config, 'memory_per_node_gb') else 256.0,
                gpus_per_node=self.config.gpus_per_node,
                is_gpu_run=is_gpu_run
            )
            
            # Get suggestions for scaling test
            suggestions = suggester.get_scaling_suggestions(
                max_nodes=self.config.max_nodes,
                scaling_type=self.config.scaling_type
            )
            
            # Display suggestion for single node
            single_node_params = suggestions[1]
            logger.info(f"")
            logger.info(f"Optimal configuration for SINGLE NODE:")
            logger.info(f"  Grid resolution:    nx={single_node_params.nx}, ny={single_node_params.ny}, nz={single_node_params.nz}")
            logger.info(f"  Particles/cell:     px={single_node_params.npcelx}, py={single_node_params.npcely}, pz={single_node_params.npcelz}")
            logger.info(f"  Decomposition:      {single_node_params.core_x}×{single_node_params.core_y}×{single_node_params.core_z} = {single_node_params.cores_used} cores")
            logger.info(f"  Memory estimate:    {single_node_params.estimated_memory_gb:.1f} GB / {single_node_params.memory_available_gb:.1f} GB")
            logger.info(f"  Core utilization:   {single_node_params.cores_used}/{single_node_params.cores_available} ({100*single_node_params.cores_used/single_node_params.cores_available:.0f}%)")
            
            # Show scaling progression if multi-node
            if self.config.max_nodes > 1:
                logger.info(f"\nScaling progression ({self.config.scaling_type.upper()})")
                logger.info(f"{'Nodes':<8} {'Grid (nx×ny×nz)':<20} {'Decomp (cx×cy×cz)':<20} {'Cores':<10} {'Memory (GB)':<15}")
                logger.info(f"-" * 75)
                
                for num_nodes in sorted(suggestions.keys()):
                    params = suggestions[num_nodes]
                    grid_str = f"{params.nx}×{params.ny}×{params.nz}"
                    decomp_str = f"{params.core_x}×{params.core_y}×{params.core_z}"
                    logger.info(f"{num_nodes:<8} {grid_str:<20} {decomp_str:<20} {params.cores_used:<10} {params.estimated_memory_gb:<15.1f}")
            
            logger.info("="*60)
            logger.info("💡 TIP: Use these values in your run.yaml or input files for optimal performance")
            logger.info("="*60)
            
        except Exception as e:
            logger.warning(f"Could not generate parameter suggestions: {e}")
    
    def _show_partition_info(self):
        """Display detailed partition hardware information."""
        try:
            # Only show if partition is specified and not a placeholder
            if self.config.partition and self.config.partition != 'X_usr_prod':
                logger.info("\n" + "="*60)
                logger.info(f"PARTITION HARDWARE SPECIFICATIONS")
                logger.info("="*60)
                
                partition_info = get_partition_info(self.config.partition)
                
                logger.info(f"Partition Name:     {self.config.partition}")
                logger.info(f"Partition Type:     {partition_info['partition_type']}")
                logger.info(f"Total Nodes:        {partition_info['total_nodes']}")
                logger.info(f"CPUs per Node:      {partition_info['cores_per_node']} cores")
                
                if partition_info['gpus_per_node'] > 0:
                    logger.info(f"GPUs per Node:      {partition_info['gpus_per_node']}")
                else:
                    logger.info(f"GPUs per Node:      0 (CPU-only partition)")
                
                logger.info(f"Memory per Node:    {partition_info['memory_per_node_gb']:.1f} GB")
                logger.info(f"Scheduler:          {partition_info['scheduler']}")
                
                # Calculate and show totals
                total_cores = partition_info['total_nodes'] * partition_info['cores_per_node']
                total_memory = partition_info['total_nodes'] * partition_info['memory_per_node_gb']
                
                logger.info(f"\nTotal Partition Resources:")
                logger.info(f"  Total CPU Cores:  {total_cores:,}")
                
                if partition_info['gpus_per_node'] > 0:
                    total_gpus = partition_info['total_nodes'] * partition_info['gpus_per_node']
                    logger.info(f"  Total GPUs:       {total_gpus:,}")
                
                logger.info(f"  Total Memory:     {total_memory:,.1f} GB")
                logger.info("="*60)
                
        except Exception as e:
            logger.debug(f"Could not retrieve partition info: {e}")
    
    def _run_tests(self) -> bool:
        """Step 4: Run scaling tests."""
        runner = TestRunner(self.test)
        success = runner.run()
        
        # Store the run directory for report generation (even if jobs not complete yet)
        self.run_dir = runner.run_dir
        self.jobs_completed = runner.jobs_completed if hasattr(runner, 'jobs_completed') else False
        
        if success and self.jobs_completed:
            logger.info("All tests completed successfully")
        elif success and not self.jobs_completed:
            logger.info("Jobs submitted successfully - waiting for completion")
        else:
            logger.error("Some tests failed or did not complete")
        
        return success
    
    def _generate_reports(self):
        """Step 5: Generate efficiency reports."""
        if not hasattr(self, 'run_dir'):
            logger.warning("No run directory available, skipping report generation")
            return
        
        try:
            generator = ReportGenerator(self.run_dir)
            report_path = generator.generate_scaling_report(
                scaling_type=self.config.scaling_type
            )
            logger.info(f"Generated scaling report: {report_path}")
            
            # Print report location
            logger.info("\n" + "="*60)
            logger.info("REPORTS GENERATED")
            logger.info("="*60)
            logger.info(f"  Text report: {report_path}")
            json_report = self.run_dir / f"{self.config.scaling_type}_scaling_report.json"
            if json_report.exists():
                logger.info(f"  JSON report: {json_report}")
            logger.info("="*60)
            
        except Exception as e:
            logger.error(f"Failed to generate reports: {e}")
    
    def _cleanup(self):
        """Clean up cloned repositories."""
        if self.source_dir and self.is_cloned:
            acquisition = CodeAcquisition(self.config.workspace_dir)
            acquisition.cleanup(self.source_dir, self.is_cloned)


def create_simple_workflow(
    source: str,
    scaling_type: str = "strong",
    max_nodes: int = 4,
    hardware_type: str = "cpu",
    **kwargs
) -> HPCOrchestrator:
    """
    Create a simple workflow with minimal configuration.
    
    Args:
        source: Local path or Git URL
        scaling_type: "strong" or "weak"
        max_nodes: Maximum number of nodes to test
        hardware_type: "cpu" or "gpu"
        **kwargs: Additional configuration options
        
    Returns:
        Configured HPCOrchestrator instance
    """
    config = OrchestratorConfig(
        source=source,
        scaling_type=scaling_type,
        max_nodes=max_nodes,
        hardware_type=hardware_type,
        **kwargs
    )
    
    return HPCOrchestrator(config)
