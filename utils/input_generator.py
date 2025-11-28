"""
Input file generator for HPC scaling tests.

This module generates input files with problem sizes and processor decompositions
tailored for specific scaling configurations.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Tuple, Any
import re

logger = logging.getLogger(__name__)


class InputFileGenerator:
    """Generate input files with scaled parameters for HPC applications."""
    
    def __init__(self, template_file: Optional[Path] = None):
        """
        Initialize input file generator.
        
        Args:
            template_file: Path to template input file (optional)
        """
        self.template_file = template_file
        self.template_content = None
        
        if template_file and template_file.exists():
            self.template_content = template_file.read_text()
            logger.info(f"Loaded template from {template_file}")
    
    def generate(self, output_file: Path, parameters: Dict[str, Any]) -> bool:
        """
        Generate input file with specified parameters.
        
        Args:
            output_file: Path where generated input file will be written
            parameters: Dictionary of parameters to substitute
            
        Returns:
            True if generation successful
        """
        if self.template_content:
            content = self._substitute_parameters(self.template_content, parameters)
        else:
            content = self._create_default_input(parameters)
        
        try:
            output_file.parent.mkdir(parents=True, exist_ok=True)
            output_file.write_text(content)
            logger.info(f"Generated input file: {output_file}")
            return True
        except Exception as e:
            logger.error(f"Failed to write input file {output_file}: {e}")
            return False
    
    def _substitute_parameters(self, content: str, parameters: Dict[str, Any]) -> str:
        """
        Substitute parameters in template content.
        
        Args:
            content: Template content
            parameters: Parameters to substitute
            
        Returns:
            Content with substituted parameters
        """
        result = content
        
        # Common parameter patterns
        param_patterns = {
            # Domain dimensions
            'Lx': r'(Lx\s*[=:]\s*)[\d.]+',
            'Ly': r'(Ly\s*[=:]\s*)[\d.]+',
            'Lz': r'(Lz\s*[=:]\s*)[\d.]+',
            
            # Grid resolution
            'nx': r'(nx\s*[=:]\s*)\d+',
            'ny': r'(ny\s*[=:]\s*)\d+',
            'nz': r'(nz\s*[=:]\s*)\d+',
            
            # Processor decomposition
            'XLEN': r'(XLEN\s*[=:]\s*)\d+',
            'YLEN': r'(YLEN\s*[=:]\s*)\d+',
            'ZLEN': r'(ZLEN\s*[=:]\s*)\d+',
            'coreX': r'(coreX\s*[=:]\s*)\d+',
            'coreY': r'(coreY\s*[=:]\s*)\d+',
            'coreZ': r'(coreZ\s*[=:]\s*)\d+',
            
            # Particles
            'npcelx': r'(npcelx\s*[=:]\s*)\d+',
            'npcely': r'(npcely\s*[=:]\s*)\d+',
            'npcelz': r'(npcelz\s*[=:]\s*)\d+',
        }
        
        for param_name, pattern in param_patterns.items():
            if param_name in parameters:
                value = parameters[param_name]
                replacement = rf'\g<1>{value}'
                result = re.sub(pattern, replacement, result, flags=re.IGNORECASE)
                logger.debug(f"Substituted {param_name} = {value}")
        
        return result
    
    def _create_default_input(self, parameters: Dict[str, Any]) -> str:
        """
        Create default input file content.
        
        Args:
            parameters: Parameters for input file
            
        Returns:
            Input file content
        """
        lines = ["# Auto-generated input file for HPC scaling test\n"]
        
        # Extract common parameters
        if 'domain_size' in parameters and parameters['domain_size']:
            dx, dy, dz = parameters['domain_size']
            lines.append(f"Lx = {dx}")
            lines.append(f"Ly = {dy}")
            lines.append(f"Lz = {dz}\n")
        
        if 'cell_count' in parameters and parameters['cell_count']:
            nx, ny, nz = parameters['cell_count']
            lines.append(f"nx = {nx}")
            lines.append(f"ny = {ny}")
            lines.append(f"nz = {nz}\n")
        
        if 'procs_decomposition' in parameters:
            px, py, pz = parameters['procs_decomposition']
            lines.append(f"# Processor decomposition")
            lines.append(f"XLEN = {px}")
            lines.append(f"YLEN = {py}")
            lines.append(f"ZLEN = {pz}\n")
        
        return '\n'.join(lines)
    
    def update_for_job_config(self, base_input: Path, output_file: Path,
                              domain_size: Optional[Tuple[float, float, float]] = None,
                              cell_count: Optional[Tuple[int, int, int]] = None,
                              procs_decomp: Optional[Tuple[int, int, int]] = None) -> bool:
        """
        Update input file for specific job configuration.
        
        Args:
            base_input: Path to base input file
            output_file: Path for output file
            domain_size: Domain dimensions (Lx, Ly, Lz)
            cell_count: Grid resolution (nx, ny, nz)
            procs_decomp: Processor decomposition (px, py, pz)
            
        Returns:
            True if successful
        """
        parameters = {}
        
        if domain_size:
            parameters['Lx'] = domain_size[0]
            parameters['Ly'] = domain_size[1]
            parameters['Lz'] = domain_size[2]
        
        if cell_count:
            parameters['nx'] = cell_count[0]
            parameters['ny'] = cell_count[1]
            parameters['nz'] = cell_count[2]
        
        if procs_decomp:
            parameters['XLEN'] = procs_decomp[0]
            parameters['YLEN'] = procs_decomp[1]
            parameters['ZLEN'] = procs_decomp[2]
            parameters['coreX'] = procs_decomp[0]
            parameters['coreY'] = procs_decomp[1]
            parameters['coreZ'] = procs_decomp[2]
        
        # Load base input if exists
        if base_input.exists():
            self.template_content = base_input.read_text()
        
        return self.generate(output_file, parameters)


def generate_input_for_scaling(base_input: Path, output_dir: Path,
                                job_configs: list) -> Dict[str, Path]:
    """
    Generate input files for all job configurations.
    
    Args:
        base_input: Base input file to use as template
        output_dir: Directory where input files will be generated
        job_configs: List of JobConfig objects
        
    Returns:
        Dictionary mapping job_id to generated input file path
    """
    generator = InputFileGenerator(base_input if base_input.exists() else None)
    generated_files = {}
    
    for job_config in job_configs:
        output_file = output_dir / f"input_{job_config.job_id}.txt"
        
        success = generator.update_for_job_config(
            base_input=base_input,
            output_file=output_file,
            domain_size=job_config.domain_size,
            cell_count=job_config.cell_count,
            procs_decomp=job_config.procs_decomposition
        )
        
        if success:
            generated_files[job_config.job_id] = output_file
    
    logger.info(f"Generated {len(generated_files)} input files in {output_dir}")
    return generated_files
