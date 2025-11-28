"""
Generic input file parser and modifier.

Completely application-agnostic parser that works with any input file format
based on user-provided variable mappings in YAML configuration.
"""

import logging
import re
from pathlib import Path
from typing import Dict, Any, Optional, Tuple

logger = logging.getLogger(__name__)


class GenericInputParser:
    """
    Parse and modify input files dynamically based on variable_map.
    
    NO hardcoded variable names - everything driven by YAML configuration.
    """
    
    def __init__(self, variable_map: Optional[Dict[str, Dict[str, str]]] = None):
        """
        Initialize generic parser.
        
        Args:
            variable_map: Mapping of parameter types to variable names
                Example:
                {
                    "length": {"x": "Lx", "y": "Ly", "z": "Lz"},
                    "cells": {"x": "nxc", "y": "nyc", "z": "nzc"},
                    "mpi": {"x": "XLEN", "y": "YLEN", "z": "ZLEN"},
                    "particles_per_cell": {"x": "npcelx", "y": "npcely", "z": "npcelz"}
                }
        """
        self.variable_map = variable_map or {}
    
    def parse_input_file(self, input_file: Path) -> Dict[str, Any]:
        """
        Parse input file and extract values for mapped variables.
        
        Args:
            input_file: Path to input file
            
        Returns:
            Dictionary of extracted values
        """
        if not input_file.exists():
            logger.warning(f"Input file not found: {input_file}")
            return {}
        
        content = input_file.read_text()
        extracted = {}
        
        # Build flat map of all variables to search for
        var_to_category = {}
        for category, dims in self.variable_map.items():
            for dim, var_name in dims.items():
                var_to_category[var_name] = (category, dim)
        
        # Search for each variable in the input file
        for var_name, (category, dim) in var_to_category.items():
            value = self._extract_variable(content, var_name)
            if value is not None:
                if category not in extracted:
                    extracted[category] = {}
                extracted[category][dim] = value
                logger.debug(f"Extracted {var_name} = {value} ({category}.{dim})")
        
        return extracted
    
    def _extract_variable(self, content: str, var_name: str) -> Optional[Any]:
        """
        Extract value of a variable from content.
        
        Supports multiple formats:
        - var_name = value
        - var_name: value
        - var_name value
        - "var_name": value
        - <var_name> value
        
        Returns numeric values as int/float, non-numeric values as strings.
        
        Args:
            content: File content
            var_name: Variable name to search for
            
        Returns:
            Extracted value (int, float, str) or None
        """
        # Build flexible regex pattern - captures ANY value (not just numeric)
        patterns = [
            rf'{re.escape(var_name)}\s*[=:]\s*([^\s#;]+)',  # var = value or var: value
            rf'["\']?{re.escape(var_name)}["\']?\s*[=:]\s*([^\s#;]+)',  # "var" = value
            rf'<{re.escape(var_name)}>\s+([^\s#;]+)',  # <var> value
            rf'{re.escape(var_name)}\s+(\S+)',  # var value (space-separated)
        ]
        
        for pattern in patterns:
            match = re.search(pattern, content, re.IGNORECASE | re.MULTILINE)
            if match:
                value_str = match.group(1).strip()
                
                # Try to convert to numeric type
                return self._safe_convert_to_numeric(value_str)
        
        return None
    
    def _safe_convert_to_numeric(self, value_str: str) -> Any:
        """
        Safely convert string to numeric type (int or float).
        
        If conversion fails, returns the original string.
        This allows placeholders like "AUTO", "PLACEHOLDER", empty strings, etc.
        
        Args:
            value_str: String value to convert
            
        Returns:
            int, float, or original string
        """
        # Handle empty or whitespace-only strings
        if not value_str or value_str.isspace():
            return value_str
        
        # Try to parse as numeric
        try:
            # Check for scientific notation or decimal point
            if '.' in value_str or 'e' in value_str.lower():
                return float(value_str)
            else:
                return int(value_str)
        except ValueError:
            # Not numeric - return as string (could be placeholder)
            return value_str
    
    def modify_input_file(
        self,
        input_file: Path,
        output_file: Path,
        modifications: Dict[str, Dict[str, Any]],
        force_replace: bool = False
    ) -> bool:
        """
        Modify input file with new parameter values.
        
        Only modifies numeric values - skips non-numeric placeholders with warning.
        
        Args:
            input_file: Source input file
            output_file: Destination for modified file
            modifications: New values in format:
                {
                    "length": {"x": 140.0, "y": 80.0},
                    "cells": {"x": 1400, "y": 800},
                    "mpi": {"x": 28, "y": 16}
                }
        
        Returns:
            True if successful
        """
        if not input_file.exists():
            logger.error(f"Input file not found: {input_file}")
            return False
        
        content = input_file.read_text()
        
        # First, extract current values to check if they're numeric
        current_values = self.parse_input_file(input_file)
        
        # Apply modifications
        for category, dims in modifications.items():
            if category not in self.variable_map:
                logger.warning(f"Unknown category '{category}' in modifications")
                continue
            
            for dim, new_value in dims.items():
                if dim not in self.variable_map[category]:
                    logger.warning(f"Unknown dimension '{dim}' in category '{category}'")
                    continue
                
                var_name = self.variable_map[category][dim]
                
                # Check if current value is numeric
                current_value = None
                if category in current_values and dim in current_values[category]:
                    current_value = current_values[category][dim]
                
                # Replace based on force_replace flag and current value type
                should_replace = False
                replace_reason = ""
                
                if current_value is None:
                    # Variable not found in file - try to replace anyway
                    should_replace = True
                    replace_reason = "variable not found"
                
                elif isinstance(current_value, (int, float)):
                    # Numeric value - safe to replace
                    should_replace = True
                    replace_reason = f"{current_value} → {new_value}"
                
                elif force_replace:
                    # Force replacement of non-numeric values (placeholders)
                    should_replace = True
                    replace_reason = f"forced replacement of '{current_value}' → {new_value}"
                
                else:
                    # Non-numeric value (placeholder/string) - skip with warning
                    logger.warning(f"  ⊗ Skipping {var_name} = '{current_value}' (non-numeric placeholder)")
                
                if should_replace:
                    old_content = content
                    content = self._replace_variable(content, var_name, new_value)
                    if content != old_content:
                        logger.info(f"  ✓ Replaced {var_name}: {replace_reason}")
                    else:
                        logger.warning(f"  ⚠ Variable {var_name} not found in input file")
        
        # Write output
        try:
            output_file.parent.mkdir(parents=True, exist_ok=True)
            output_file.write_text(content)
            logger.info(f"Generated modified input: {output_file}")
            return True
        except Exception as e:
            logger.error(f"Failed to write output file: {e}")
            return False
    
    def _replace_variable(self, content: str, var_name: str, new_value: Any) -> str:
        """
        Replace variable value in content.
        
        Handles both numeric and non-numeric values in the input file.
        
        Args:
            content: File content
            var_name: Variable name
            new_value: New value
            
        Returns:
            Modified content
        """
        # Build replacement patterns - match ANY value (not just numeric)
        # Use word boundary or non-word character to avoid partial matches
        patterns = [
            # Pattern 1: var = value or var: value
            (rf'({re.escape(var_name)}\s*[=:])\s*\S+', rf'\g<1> {new_value}'),
            # Pattern 2: "var" = value or 'var': value  
            (rf'(["\']?{re.escape(var_name)}["\']?\s*[=:])\s*\S+', rf'\g<1> {new_value}'),
            # Pattern 3: <var> value
            (rf'(<{re.escape(var_name)}>)\s+\S+', rf'\g<1> {new_value}'),
        ]
        
        for pattern, replacement in patterns:
            new_content = re.sub(pattern, replacement, content, flags=re.IGNORECASE | re.MULTILINE, count=1)
            if new_content != content:
                return new_content
        
        # If no pattern matched, log warning but don't fail
        logger.debug(f"Could not find pattern to replace '{var_name}' in input file")
        return content
    
    def generate_scaled_input(
        self,
        base_input: Path,
        output_file: Path,
        domain_size: Optional[Tuple[float, float, float]] = None,
        cell_count: Optional[Tuple[int, int, int]] = None,
        mpi_decomp: Optional[Tuple[int, int, int]] = None,
        particles_per_cell: Optional[Tuple[int, int, int]] = None
    ) -> bool:
        """
        Generate scaled input file for a specific configuration.
        
        Args:
            base_input: Base input file
            output_file: Output path
            domain_size: (Lx, Ly, Lz) values
            cell_count: (nx, ny, nz) values
            mpi_decomp: (px, py, pz) values
            particles_per_cell: (npcelx, npcely, npcelz) values
            
        Returns:
            True if successful
        """
        logger.info(f"═══ GenericInputParser.generate_scaled_input ===")
        logger.info(f"Base input: {base_input}")
        logger.info(f"Output file: {output_file}")
        logger.info(f"Domain size: {domain_size}")
        logger.info(f"Cell count: {cell_count}")
        logger.info(f"MPI decomp: {mpi_decomp}")
        logger.info(f"Particles/cell: {particles_per_cell}")
        logger.info(f"Variable map: {self.variable_map}")
        logger.info(f"═══════════════════════════════════════════════")
        
        modifications = {}
        
        if domain_size and "length" in self.variable_map:
            modifications["length"] = {
                "x": domain_size[0],
                "y": domain_size[1],
                "z": domain_size[2]
            }
        
        if cell_count and "cells" in self.variable_map:
            modifications["cells"] = {
                "x": cell_count[0],
                "y": cell_count[1],
                "z": cell_count[2]
            }
        
        if mpi_decomp and ("mpi" in self.variable_map or "mpi_decomposition" in self.variable_map):
            # Support both 'mpi' and 'mpi_decomposition' as category names
            mpi_category = "mpi" if "mpi" in self.variable_map else "mpi_decomposition"
            modifications[mpi_category] = {
                "x": mpi_decomp[0],
                "y": mpi_decomp[1],
                "z": mpi_decomp[2]
            }
        
        if particles_per_cell and ("particles_per_cell" in self.variable_map or "particles" in self.variable_map):
            # Support both 'particles_per_cell' and 'particles' as category names
            particles_category = "particles_per_cell" if "particles_per_cell" in self.variable_map else "particles"
            modifications[particles_category] = {
                "x": particles_per_cell[0],
                "y": particles_per_cell[1],
                "z": particles_per_cell[2]
            }
        
        return self.modify_input_file(base_input, output_file, modifications, force_replace=True)
