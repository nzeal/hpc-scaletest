"""Intelligent Input File Analyzer and Validator

Parses, validates, and adapts simulation input parameters for scaling tests.
Ensures consistency between grid resolution, domain decomposition, and scaling configuration.
Supports flexible variable name mapping for different simulation codes.
"""

import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


# Default parameter mapping for common simulation codes
DEFAULT_PARAMETER_MAP = {
    # Domain size (physical dimensions)
    'domain_x': ['Lx', 'LX', 'DomainX', 'XLength', 'LENGTH_X', 'domain_size_x', 'L_x'],
    'domain_y': ['Ly', 'LY', 'DomainY', 'YLength', 'LENGTH_Y', 'domain_size_y', 'L_y'],
    'domain_z': ['Lz', 'LZ', 'DomainZ', 'ZLength', 'LENGTH_Z', 'domain_size_z', 'L_z'],
    
    # Grid resolution (cell counts)
    'grid_x': ['nxc', 'nx', 'NX', 'GridX', 'CellsX', 'NCELLS_X', 'grid_size_x', 'n_x'],
    'grid_y': ['nyc', 'ny', 'NY', 'GridY', 'CellsY', 'NCELLS_Y', 'grid_size_y', 'n_y'],
    'grid_z': ['nzc', 'nz', 'NZ', 'GridZ', 'CellsZ', 'NCELLS_Z', 'grid_size_z', 'n_z'],
    
    # MPI decomposition (processor counts)
    'decomp_x': ['XLEN', 'nprocx', 'NPROCX', 'DecompX', 'coreX', 'nproc_x', 'px'],
    'decomp_y': ['YLEN', 'nprocy', 'NPROCY', 'DecompY', 'coreY', 'nproc_y', 'py'],
    'decomp_z': ['ZLEN', 'nprocz', 'NPROCZ', 'DecompZ', 'coreZ', 'nproc_z', 'pz'],
    
    # Particles per cell
    'particles_x': ['npcelx', 'ParticlesX', 'npartx', 'Npx', 'particles_per_cell_x'],
    'particles_y': ['npcely', 'ParticlesY', 'nparty', 'Npy', 'particles_per_cell_y'],
    'particles_z': ['npcelz', 'ParticlesZ', 'npartz', 'Npz', 'particles_per_cell_z'],
    
    # Time stepping
    'timestep': ['dt', 'DT', 'TimeStep', 'time_step', 'delta_t'],
    'num_cycles': ['ncycles', 'NCYCLES', 'NumCycles', 'num_steps', 'nsteps'],
}


@dataclass
class SimulationParameters:
    """Parsed simulation parameters from input file."""
    # Domain size
    Lx: Optional[float] = None
    Ly: Optional[float] = None
    Lz: Optional[float] = None
    
    # Grid resolution
    nxc: Optional[int] = None
    nyc: Optional[int] = None
    nzc: Optional[int] = None
    
    # Domain decomposition
    XLEN: Optional[int] = None
    YLEN: Optional[int] = None
    ZLEN: Optional[int] = None
    
    # Particles per cell
    npcelx: Optional[int] = None
    npcely: Optional[int] = None
    npcelz: Optional[int] = None
    
    # Species-specific particles
    npcel: Dict[str, int] = field(default_factory=dict)
    
    # Other parameters
    dt: Optional[float] = None
    ncycles: Optional[int] = None
    
    # Raw content
    raw_lines: List[str] = field(default_factory=list)
    parameter_map: Dict[str, int] = field(default_factory=dict)  # param_name -> line_index


@dataclass
class ValidationResult:
    """Result of parameter validation."""
    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    suggestions: List[str] = field(default_factory=list)


class InputFileAnalyzer:
    """
    Intelligent input file analyzer with parameter validation and adaptation.
    
    Capabilities:
    - Parse simulation input files
    - Validate parameter consistency
    - Detect decomposition/grid mismatches
    - Suggest optimal configurations
    - Dynamically adapt parameters for scaling
    """
    
    def __init__(self, input_file: Path, parameter_map: Optional[Dict[str, List[str]]] = None):
        """
        Initialize analyzer with input file.
        
        Args:
            input_file: Path to simulation input file
            parameter_map: Custom parameter name mapping (uses DEFAULT_PARAMETER_MAP if None)
        """
        self.input_file = input_file
        self.params = SimulationParameters()
        self.validation = ValidationResult(is_valid=True)
        
        # Use custom mapping or default
        self.param_map = parameter_map or DEFAULT_PARAMETER_MAP
        
        # Build reverse lookup: variable_name -> canonical_name
        self.name_to_canonical = {}
        for canonical, variants in self.param_map.items():
            for variant in variants:
                self.name_to_canonical[variant] = canonical
        
        if input_file and input_file.exists():
            self._parse_input_file()
    
    def _parse_input_file(self):
        """Parse input file and extract parameters."""
        logger.info(f"Parsing input file: {self.input_file}")
        
        try:
            with open(self.input_file, 'r') as f:
                self.params.raw_lines = f.readlines()
            
            for idx, line in enumerate(self.params.raw_lines):
                stripped = line.strip()
                
                # Skip comments and empty lines
                if not stripped or stripped.startswith('#') or stripped.startswith('//'):
                    continue
                
                # Parse parameter = value patterns
                if '=' in stripped:
                    parts = stripped.split('=', 1)
                    if len(parts) == 2:
                        param_name = parts[0].strip()
                        param_value = parts[1].strip().split()[0]  # Get first token
                        
                        self._extract_parameter(param_name, param_value, idx)
            
            logger.info(f"✓ Parsed {len([p for p in self.params.parameter_map])} parameters")
            self._log_parsed_params()
            
        except Exception as e:
            logger.error(f"Failed to parse input file: {e}")
            raise
    
    def _extract_parameter(self, name: str, value: str, line_idx: int):
        """Extract and store parameter value using flexible name mapping."""
        self.params.parameter_map[name] = line_idx
        
        # Try to find canonical name from mapping
        canonical = self.name_to_canonical.get(name)
        
        if not canonical:
            # Not in mapping - ignore or log
            logger.debug(f"Skipping unknown parameter: {name}")
            return
        
        try:
            # Map to internal parameter based on canonical name
            if canonical == 'domain_x':
                self.params.Lx = float(value)
            elif canonical == 'domain_y':
                self.params.Ly = float(value)
            elif canonical == 'domain_z':
                self.params.Lz = float(value)
            
            elif canonical == 'grid_x':
                self.params.nxc = int(value)
            elif canonical == 'grid_y':
                self.params.nyc = int(value)
            elif canonical == 'grid_z':
                self.params.nzc = int(value)
            
            elif canonical == 'decomp_x':
                self.params.XLEN = int(value)
            elif canonical == 'decomp_y':
                self.params.YLEN = int(value)
            elif canonical == 'decomp_z':
                self.params.ZLEN = int(value)
            
            elif canonical == 'particles_x':
                self.params.npcelx = int(value)
            elif canonical == 'particles_y':
                self.params.npcely = int(value)
            elif canonical == 'particles_z':
                self.params.npcelz = int(value)
            
            elif canonical == 'timestep':
                self.params.dt = float(value)
            elif canonical == 'num_cycles':
                self.params.ncycles = int(value)
                
            logger.debug(f"  Mapped {name} -> {canonical} = {value}")
            
        except ValueError as e:
            logger.warning(f"Could not parse {name}={value}: {e}")
    
    def _log_parsed_params(self):
        """Log parsed parameters for debugging."""
        logger.debug("Parsed parameters:")
        if self.params.Lx is not None:
            logger.debug(f"  Domain: Lx={self.params.Lx}, Ly={self.params.Ly}, Lz={self.params.Lz}")
        if self.params.nxc is not None:
            logger.debug(f"  Grid: nxc={self.params.nxc}, nyc={self.params.nyc}, nzc={self.params.nzc}")
        if self.params.XLEN is not None:
            logger.debug(f"  Decomp: XLEN={self.params.XLEN}, YLEN={self.params.YLEN}, ZLEN={self.params.ZLEN}")
    
    def validate(self, num_procs: int) -> ValidationResult:
        """
        Validate parameter consistency.
        
        Args:
            num_procs: Total number of MPI processes
            
        Returns:
            ValidationResult with errors, warnings, and suggestions
        """
        self.validation = ValidationResult(is_valid=True)
        
        logger.info(f"\n{'='*60}")
        logger.info("VALIDATING INPUT PARAMETERS")
        logger.info(f"{'='*60}")
        
        # Check decomposition vs total processes
        self._validate_decomposition(num_procs)
        
        # Check grid divisibility by decomposition
        self._validate_grid_decomposition()
        
        # Check domain/grid consistency
        self._validate_domain_grid_consistency()
        
        # Check physical constraints
        self._validate_physical_constraints()
        
        # Generate report
        self._generate_validation_report()
        
        return self.validation
    
    def _validate_decomposition(self, num_procs: int):
        """Validate decomposition matches process count."""
        if self.params.XLEN and self.params.YLEN and self.params.ZLEN:
            expected_procs = self.params.XLEN * self.params.YLEN * self.params.ZLEN
            
            if expected_procs != num_procs:
                error = f"Decomposition mismatch: XLEN({self.params.XLEN}) × YLEN({self.params.YLEN}) × ZLEN({self.params.ZLEN}) = {expected_procs}, but num_procs = {num_procs}"
                self.validation.errors.append(error)
                self.validation.is_valid = False
                
                # Suggest optimal decomposition
                suggestion = self._suggest_optimal_decomposition(num_procs)
                self.validation.suggestions.append(f"Suggested decomposition: {suggestion}")
        else:
            self.validation.warnings.append("Decomposition parameters (XLEN, YLEN, ZLEN) not found in input file")
    
    def _validate_grid_decomposition(self):
        """Validate grid is divisible by decomposition."""
        if self.params.nxc and self.params.XLEN:
            if self.params.nxc % self.params.XLEN != 0:
                error = f"Grid X not divisible: nxc({self.params.nxc}) % XLEN({self.params.XLEN}) = {self.params.nxc % self.params.XLEN}"
                self.validation.errors.append(error)
                self.validation.is_valid = False
                self.validation.suggestions.append(f"Adjust nxc to multiple of {self.params.XLEN}, e.g., {(self.params.nxc // self.params.XLEN) * self.params.XLEN}")
        
        if self.params.nyc and self.params.YLEN:
            if self.params.nyc % self.params.YLEN != 0:
                error = f"Grid Y not divisible: nyc({self.params.nyc}) % YLEN({self.params.YLEN}) = {self.params.nyc % self.params.YLEN}"
                self.validation.errors.append(error)
                self.validation.is_valid = False
                self.validation.suggestions.append(f"Adjust nyc to multiple of {self.params.YLEN}, e.g., {(self.params.nyc // self.params.YLEN) * self.params.YLEN}")
        
        if self.params.nzc and self.params.ZLEN:
            if self.params.nzc % self.params.ZLEN != 0:
                error = f"Grid Z not divisible: nzc({self.params.nzc}) % ZLEN({self.params.ZLEN}) = {self.params.nzc % self.params.ZLEN}"
                self.validation.errors.append(error)
                self.validation.is_valid = False
                self.validation.suggestions.append(f"Adjust nzc to multiple of {self.params.ZLEN}, e.g., {(self.params.nzc // self.params.ZLEN + 1) * self.params.ZLEN}")
    
    def _validate_domain_grid_consistency(self):
        """Validate domain size and grid resolution are consistent."""
        if self.params.Lx and self.params.nxc:
            dx = self.params.Lx / self.params.nxc
            if self.params.Ly and self.params.nyc:
                dy = self.params.Ly / self.params.nyc
                if abs(dx - dy) > 0.1 * max(dx, dy):  # More than 10% difference
                    self.validation.warnings.append(f"Non-uniform grid spacing: dx={dx:.4f}, dy={dy:.4f}")
    
    def _validate_physical_constraints(self):
        """Validate physical constraints."""
        # Check for very small grid cells
        if self.params.Lx and self.params.nxc:
            cell_size_x = self.params.Lx / self.params.nxc
            if cell_size_x < 0.001:
                self.validation.warnings.append(f"Very small cell size in X: {cell_size_x}")
        
        # Check for reasonable timestep
        if self.params.dt:
            if self.params.dt > 1.0:
                self.validation.warnings.append(f"Large timestep: dt={self.params.dt} (may cause instability)")
            elif self.params.dt < 1e-6:
                self.validation.warnings.append(f"Very small timestep: dt={self.params.dt} (may be slow)")
    
    def _suggest_optimal_decomposition(self, num_procs: int) -> str:
        """
        Suggest optimal decomposition considering grid dimensions.
        
        Args:
            num_procs: Total number of processes
            
        Returns:
            Suggested decomposition string
        """
        logger.info(f"Computing optimal decomposition for {num_procs} processes")
        
        # Use intelligent grid-aware decomposition
        best_x, best_y, best_z = self._compute_grid_aware_decomposition(num_procs)
        
        return f"XLEN={best_x}, YLEN={best_y}, ZLEN={best_z}"
    
    def _compute_grid_aware_decomposition(self, num_procs: int) -> Tuple[int, int, int]:
        """
        Compute optimal 3D decomposition considering grid dimensions.
        
        Strategy:
        1. Prioritize splitting along directions with larger grid counts
        2. Ensure no dimension gets more splits than grid cells
        3. Maintain good aspect ratio for communication efficiency
        
        Args:
            num_procs: Total number of processes
            
        Returns:
            Optimal decomposition (nprocx, nprocy, nprocz)
        """
        # Get grid dimensions
        nx = self.params.nxc if self.params.nxc else 100
        ny = self.params.nyc if self.params.nyc else 100
        nz = self.params.nzc if self.params.nzc else 1
        
        logger.info(f"Grid dimensions: nx={nx}, ny={ny}, nz={nz}")
        
        # Sort dimensions by grid size (largest first)
        dims = [
            ('x', nx, self.params.Lx if self.params.Lx else nx),
            ('y', ny, self.params.Ly if self.params.Ly else ny),
            ('z', nz, self.params.Lz if self.params.Lz else nz)
        ]
        dims.sort(key=lambda d: d[1], reverse=True)  # Sort by grid count
        
        logger.info(f"Dimension priority: {[f'{d[0]}({d[1]} cells)' for d in dims]}")
        
        # Find all valid factorizations
        valid_decomps = []
        
        for px in range(1, min(nx, num_procs) + 1):
            if num_procs % px != 0:
                continue
            remaining = num_procs // px
            
            for py in range(1, min(ny, remaining) + 1):
                if remaining % py != 0:
                    continue
                pz = remaining // py
                
                # Check constraints
                if pz > nz:  # Can't split more than grid cells
                    continue
                
                if px * py * pz != num_procs:
                    continue
                
                # Score this decomposition
                score = self._score_decomposition(px, py, pz, nx, ny, nz)
                valid_decomps.append((score, px, py, pz))
        
        if not valid_decomps:
            logger.warning("No valid decomposition found, using fallback")
            return self._factorize_3d(num_procs, (1.0, 1.0, 1.0))
        
        # Choose best decomposition
        valid_decomps.sort(key=lambda x: x[0])
        best_score, best_px, best_py, best_pz = valid_decomps[0]
        
        logger.info(f"Selected decomposition: {best_px}×{best_py}×{best_pz} (score={best_score:.3f})")
        logger.info(f"  Cells per subdomain: X={nx//best_px}, Y={ny//best_py}, Z={nz//best_pz}")
        
        return (best_px, best_py, best_pz)
    
    def _score_decomposition(self, px: int, py: int, pz: int, 
                            nx: int, ny: int, nz: int) -> float:
        """
        Score a decomposition based on multiple criteria.
        Lower score is better.
        
        Criteria:
        1. Balance of subdomain sizes
        2. Minimization of surface area (communication)
        3. Preference for splitting larger dimensions
        
        Args:
            px, py, pz: Processor counts in each direction
            nx, ny, nz: Grid cells in each direction
            
        Returns:
            Score (lower is better)
        """
        # Cells per subdomain in each direction
        cells_x = nx / px if px > 0 else 0
        cells_y = ny / py if py > 0 else 0
        cells_z = nz / pz if pz > 0 else 0
        
        # Avoid division by zero
        if cells_x == 0 or cells_y == 0 or cells_z == 0:
            return float('inf')
        
        # 1. Balance score - prefer similar subdomain sizes
        max_cells = max(cells_x, cells_y, cells_z)
        min_cells = min(cells_x, cells_y, cells_z)
        balance_score = (max_cells - min_cells) / max_cells if max_cells > 0 else 0
        
        # 2. Communication score - prefer minimal surface area
        # Surface area proportional to faces between subdomains
        surface_x = py * pz * 2  # Faces perpendicular to X
        surface_y = px * pz * 2  # Faces perpendicular to Y
        surface_z = px * py * 2  # Faces perpendicular to Z
        total_surface = surface_x + surface_y + surface_z
        
        volume = px * py * pz
        surface_to_volume = total_surface / volume if volume > 0 else float('inf')
        
        # 3. Anisotropy penalty - avoid extreme ratios
        max_procs = max(px, py, pz)
        min_procs = min(px, py, pz)
        anisotropy = (max_procs / min_procs) if min_procs > 0 else float('inf')
        
        # 4. Grid fit penalty - penalize splitting dimensions with few cells
        grid_fit_penalty = 0
        if nz == 1 and pz > 1:
            grid_fit_penalty += 1000  # Severe penalty for splitting 1-cell dimension
        if cells_z < 2 and pz > 1:
            grid_fit_penalty += 100   # High penalty for very small subdomains
        if cells_x < 2 or cells_y < 2:
            grid_fit_penalty += 50
        
        # Combined score (weighted)
        score = (balance_score * 2.0 +           # Balance is important
                surface_to_volume * 0.5 +        # Communication overhead
                (anisotropy - 1.0) * 0.3 +       # Avoid extreme ratios
                grid_fit_penalty)                 # Respect grid constraints
        
        return score
    
    def _factorize_3d(self, n: int, aspect_ratio: Tuple[float, float, float]) -> Tuple[int, int, int]:
        """Find optimal 3D factorization maintaining aspect ratio."""
        best_score = float('inf')
        best_decomp = (1, 1, n)
        
        for x in range(1, int(n**0.5) + 1):
            if n % x != 0:
                continue
            remaining = n // x
            for y in range(1, int(remaining**0.5) + 1):
                if remaining % y != 0:
                    continue
                z = remaining // y
                
                # Score based on deviation from aspect ratio
                total = x + y + z
                ratio_x = x / total if total > 0 else 0
                ratio_y = y / total if total > 0 else 0
                ratio_z = z / total if total > 0 else 0
                
                score = (abs(ratio_x - aspect_ratio[0]) +
                        abs(ratio_y - aspect_ratio[1]) +
                        abs(ratio_z - aspect_ratio[2]))
                
                if score < best_score:
                    best_score = score
                    best_decomp = (x, y, z)
        
        return best_decomp
    
    def _generate_validation_report(self):
        """Generate validation report."""
        if self.validation.errors:
            logger.error(f"\n❌ VALIDATION FAILED - {len(self.validation.errors)} error(s) found:")
            for i, error in enumerate(self.validation.errors, 1):
                logger.error(f"  {i}. {error}")
        
        if self.validation.warnings:
            logger.warning(f"\n⚠ {len(self.validation.warnings)} warning(s):")
            for i, warning in enumerate(self.validation.warnings, 1):
                logger.warning(f"  {i}. {warning}")
        
        if self.validation.suggestions:
            logger.info(f"\n💡 Suggestions:")
            for i, suggestion in enumerate(self.validation.suggestions, 1):
                logger.info(f"  {i}. {suggestion}")
        
        if self.validation.is_valid and not self.validation.warnings:
            logger.info(f"\n✅ VALIDATION PASSED - All parameters are consistent")
        
        logger.info(f"{'='*60}\n")
    
    def adapt_for_scaling(
        self,
        scaling_type: str,
        num_procs: int,
        procs_decomp: Optional[Tuple[int, int, int]] = None,
        scale_factor: float = 1.0
    ) -> 'SimulationParameters':
        """
        Adapt parameters for scaling test with intelligent decomposition.
        
        Args:
            scaling_type: 'strong' or 'weak'
            num_procs: Total number of processes
            procs_decomp: Processor decomposition (px, py, pz), None for auto
            scale_factor: Scaling factor for weak scaling
            
        Returns:
            Adapted SimulationParameters
        """
        logger.info(f"Adapting parameters for {scaling_type} scaling (scale_factor={scale_factor})")
        
        adapted = SimulationParameters()
        adapted.raw_lines = self.params.raw_lines.copy()
        adapted.parameter_map = self.params.parameter_map.copy()
        
        # Compute optimal decomposition if not provided
        if procs_decomp is None:
            logger.info("Auto-computing optimal decomposition...")
            procs_decomp = self._compute_grid_aware_decomposition(num_procs)
        else:
            logger.info(f"Using provided decomposition: {procs_decomp}")
        
        # Update decomposition
        adapted.XLEN, adapted.YLEN, adapted.ZLEN = procs_decomp
        
        if scaling_type == 'weak':
            # Weak scaling: increase problem size proportionally
            linear_scale = scale_factor ** (1.0 / 3.0)  # Cubic root for 3D
            
            # Scale grid - ensure divisibility by decomposition
            if self.params.nxc:
                scaled_nx = int(self.params.nxc * linear_scale)
                # Round to nearest multiple of XLEN
                adapted.nxc = max(adapted.XLEN, (scaled_nx // adapted.XLEN) * adapted.XLEN)
            
            if self.params.nyc:
                scaled_ny = int(self.params.nyc * linear_scale)
                adapted.nyc = max(adapted.YLEN, (scaled_ny // adapted.YLEN) * adapted.YLEN)
            
            if self.params.nzc:
                scaled_nz = int(self.params.nzc * linear_scale)
                adapted.nzc = max(adapted.ZLEN, (scaled_nz // adapted.ZLEN) * adapted.ZLEN)
            
            # Scale domain proportionally
            if self.params.Lx and adapted.nxc and self.params.nxc:
                adapted.Lx = self.params.Lx * (adapted.nxc / self.params.nxc)
            if self.params.Ly and adapted.nyc and self.params.nyc:
                adapted.Ly = self.params.Ly * (adapted.nyc / self.params.nyc)
            if self.params.Lz and adapted.nzc and self.params.nzc:
                adapted.Lz = self.params.Lz * (adapted.nzc / self.params.nzc)
        
        else:  # strong scaling
            # Strong scaling: keep problem size constant, adjust decomposition only
            adapted.nxc = self.params.nxc
            adapted.nyc = self.params.nyc
            adapted.nzc = self.params.nzc
            adapted.Lx = self.params.Lx
            adapted.Ly = self.params.Ly
            adapted.Lz = self.params.Lz
            
            # Ensure grid is divisible by new decomposition
            if adapted.nxc and adapted.nxc % adapted.XLEN != 0:
                logger.warning(f"Grid X ({adapted.nxc}) not divisible by XLEN ({adapted.XLEN})")
                adapted.nxc = (adapted.nxc // adapted.XLEN + 1) * adapted.XLEN
                logger.info(f"  Adjusted nxc to {adapted.nxc}")
            
            if adapted.nyc and adapted.nyc % adapted.YLEN != 0:
                logger.warning(f"Grid Y ({adapted.nyc}) not divisible by YLEN ({adapted.YLEN})")
                adapted.nyc = (adapted.nyc // adapted.YLEN + 1) * adapted.YLEN
                logger.info(f"  Adjusted nyc to {adapted.nyc}")
            
            if adapted.nzc and adapted.nzc % adapted.ZLEN != 0:
                logger.warning(f"Grid Z ({adapted.nzc}) not divisible by ZLEN ({adapted.ZLEN})")
                adapted.nzc = (adapted.nzc // adapted.ZLEN + 1) * adapted.ZLEN
                logger.info(f"  Adjusted nzc to {adapted.nzc}")
        
        # Keep particles per cell constant
        adapted.npcelx = self.params.npcelx
        adapted.npcely = self.params.npcely
        adapted.npcelz = self.params.npcelz
        
        # Keep timestep and cycles
        adapted.dt = self.params.dt
        adapted.ncycles = self.params.ncycles
        
        logger.info(f"  Decomposition: {adapted.XLEN}×{adapted.YLEN}×{adapted.ZLEN}")
        if adapted.nxc:
            logger.info(f"  Grid: {adapted.nxc}×{adapted.nyc}×{adapted.nzc}")
        if adapted.Lx:
            logger.info(f"  Domain: {adapted.Lx:.2f}×{adapted.Ly:.2f}×{adapted.Lz:.2f}")
        
        return adapted
    
    def generate_adapted_input(self, adapted_params: SimulationParameters, output_path: Path):
        """
        Generate new input file with adapted parameters.
        
        Args:
            adapted_params: Adapted simulation parameters
            output_path: Path for output file
        """
        logger.info(f"Generating adapted input file: {output_path}")
        
        new_lines = []
        for idx, line in enumerate(adapted_params.raw_lines):
            # Check if this line contains a parameter we need to update
            updated = False
            
            for param_name, param_line_idx in adapted_params.parameter_map.items():
                if idx == param_line_idx:
                    # Get new value
                    new_value = getattr(adapted_params, param_name, None)
                    if new_value is not None:
                        # Preserve formatting, just update value
                        if '=' in line:
                            parts = line.split('=', 1)
                            new_line = f"{parts[0]}= {new_value}\n"
                            new_lines.append(new_line)
                            updated = True
                            logger.debug(f"  Updated {param_name} = {new_value}")
                            break
            
            if not updated:
                new_lines.append(line)
        
        # Write to file
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            f.writelines(new_lines)
        
        logger.info(f"✓ Generated adapted input file")


def validate_and_adapt_input(
    input_file: Path,
    num_procs: int,
    procs_decomp: Optional[Tuple[int, int, int]] = None,
    scaling_type: str = 'strong',
    scale_factor: float = 1.0,
    output_file: Optional[Path] = None,
    use_llm: bool = False,
    llm_api_key: Optional[str] = None,
    llm_model: str = "gpt-4",
    custom_param_map: Optional[Dict[str, List[str]]] = None
) -> Tuple[bool, Optional[Path]]:
    """
    Validate and adapt input file for scaling test.
    Uses intelligent grid-aware decomposition if procs_decomp is None.
    Optionally uses LLM for parameter discovery.
    
    Args:
        input_file: Original input file
        num_procs: Total number of processes
        procs_decomp: Processor decomposition (px, py, pz), None for auto
        scaling_type: 'strong' or 'weak'
        scale_factor: Scaling factor for weak scaling
        output_file: Path for adapted input file
        use_llm: Enable LLM-enhanced parameter discovery
        llm_api_key: OpenAI API key (optional if env var set)
        llm_model: LLM model to use (default: gpt-4)
        custom_param_map: Custom parameter name mappings from YAML
        
    Returns:
        Tuple of (is_valid, adapted_file_path)
    """
    # Create analyzer with LLM support if enabled
    if use_llm or custom_param_map:
        try:
            from utils.llm_parameter_mapper import create_llm_enhanced_analyzer
            
            logger.info("Creating LLM-enhanced input analyzer...")
            analyzer = create_llm_enhanced_analyzer(
                input_file=input_file,
                use_llm=use_llm,
                api_key=llm_api_key,
                model=llm_model,
                yaml_aliases=custom_param_map
            )
            logger.info("✓ LLM-enhanced analyzer created")
        except ImportError as e:
            logger.warning(f"LLM mapper not available: {e}. Using rule-based analyzer.")
            analyzer = InputFileAnalyzer(input_file, parameter_map=custom_param_map)
    else:
        # Standard rule-based analyzer
        analyzer = InputFileAnalyzer(input_file, parameter_map=custom_param_map)
    
    # Validate original parameters
    validation = analyzer.validate(num_procs)
    
    if not validation.is_valid:
        logger.error("Input file validation failed. Attempting auto-correction...")
    
    # Adapt parameters
    adapted_params = analyzer.adapt_for_scaling(
        scaling_type=scaling_type,
        num_procs=num_procs,
        procs_decomp=procs_decomp,
        scale_factor=scale_factor
    )
    
    # Validate adapted parameters
    temp_analyzer = InputFileAnalyzer(input_file)  # Reuse input_file instead of None
    temp_analyzer.params = adapted_params
    adapted_validation = temp_analyzer.validate(num_procs)
    
    if not adapted_validation.is_valid:
        logger.error("Adapted parameters still invalid!")
        return False, None
    
    # Generate adapted file
    if output_file:
        analyzer.generate_adapted_input(adapted_params, output_file)
        return True, output_file
    
    return True, None
