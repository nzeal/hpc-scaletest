"""
LLM-Enhanced Parameter Mapping for Unknown Simulation Codes

This module provides OPTIONAL LLM-based parameter discovery for input files
with completely unknown variable names. It ENHANCES (not replaces) the
rule-based mapping system.

Usage:
    # Automatic fallback to LLM if needed
    analyzer = InputFileAnalyzer(input_file, use_llm_fallback=True)
    
    # Or explicit LLM mapping generation
    mapper = LLMParameterMapper()
    custom_map = mapper.discover_parameters(input_file)
    analyzer = InputFileAnalyzer(input_file, parameter_map=custom_map)
"""

import logging
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re

logger = logging.getLogger(__name__)

# Try to import OpenAI (optional dependency)
try:
    import openai
    HAS_OPENAI = True
except ImportError:
    HAS_OPENAI = False
    logger.debug("OpenAI not available - LLM features disabled")


class LLMParameterMapper:
    """
    LLM-powered parameter mapper for unknown simulation codes.
    
    Uses AI to discover parameter meanings from context, comments,
    and variable names in input files.
    """
    
    def __init__(self, api_key: Optional[str] = None, model: str = "gpt-4"):
        """
        Initialize LLM parameter mapper.
        
        Args:
            api_key: OpenAI API key (optional, reads from env if not provided)
            model: LLM model to use (default: gpt-4)
        """
        self.api_key = api_key
        self.model = model
        self.enabled = HAS_OPENAI
        
        if not self.enabled:
            logger.warning("LLM features disabled - install openai package to enable")
    
    def discover_parameters(
        self,
        input_file: Path,
        sample_lines: int = 100
    ) -> Dict[str, List[str]]:
        """
        Discover parameter mappings using LLM analysis.
        
        Args:
            input_file: Path to input file
            sample_lines: Number of lines to analyze (default: 100)
            
        Returns:
            Parameter mapping dict (canonical_name -> [variant_names])
        """
        if not self.enabled:
            logger.warning("LLM discovery not available - using default mapping")
            from utils.input_analyzer import DEFAULT_PARAMETER_MAP
            return DEFAULT_PARAMETER_MAP
        
        logger.info(f"Discovering parameters in {input_file.name} using LLM...")
        
        # Read sample of input file
        with open(input_file, 'r') as f:
            lines = f.readlines()[:sample_lines]
        
        sample_content = ''.join(lines)
        
        # Build prompt for LLM
        prompt = self._build_discovery_prompt(sample_content)
        
        try:
            # Call LLM API
            response = self._call_llm(prompt)
            
            # Parse response
            mapping = self._parse_llm_response(response)
            
            logger.info(f"✓ Discovered {len(mapping)} parameter types")
            return mapping
            
        except Exception as e:
            logger.error(f"LLM discovery failed: {e}")
            logger.info("Falling back to default parameter mapping")
            from utils.input_analyzer import DEFAULT_PARAMETER_MAP
            return DEFAULT_PARAMETER_MAP
    
    def _build_discovery_prompt(self, content: str) -> str:
        """Build prompt for LLM to discover parameters."""
        prompt = f"""You are an expert in HPC simulation input files. Analyze this simulation input file and identify key parameters.

INPUT FILE SAMPLE:
```
{content}
```

TASK:
Identify which variables correspond to these canonical parameter types by analyzing:
1. Variable names (exact matches, similar names, abbreviations)
2. Comments explaining what each parameter means
3. Context clues (units, typical values, grouping with related parameters)
4. Physical meanings (domain size vs grid resolution vs decomposition)

CANONICAL PARAMETER TYPES:

**Domain Size (Physical Box Dimensions):**
- domain_x, domain_y, domain_z
- Physical length of simulation box in each direction
- Common names: Lx, Ly, Lz, DomainX, BoxLength, XLength, LENGTH_X, L_x
- Units: typically meters, cm, or dimensionless
- Values: positive real numbers

**Grid Resolution (Number of Cells):**
- grid_x, grid_y, grid_z
- Number of grid cells/points in each direction
- Common names: nx, ny, nz, nxc, nyc, nzc, GridX, CellsX, NCELLS_X, n_x
- Values: positive integers (often powers of 2)

**MPI Decomposition (Processor Counts):**
- decomp_x, decomp_y, decomp_z
- Number of MPI processes/subdomains in each direction
- Common names: XLEN, YLEN, ZLEN, nprocx, nprocy, nprocz, DecompX, coreX, px, py, pz
- Values: positive integers
- Constraint: product must equal total MPI processes

**Particles Per Cell:**
- particles_x, particles_y, particles_z
- Number of particles per grid cell in each direction
- Common names: npcelx, npcely, npcelz, ParticlesX, npartx, ppc_x
- Values: positive integers

**Time Stepping:**
- timestep: Time step size (dt, DT, TimeStep, delta_t)
- num_cycles: Number of time steps (ncycles, NCYCLES, nsteps, NumCycles)

ANALYSIS STRATEGY:
1. Look for exact name matches first
2. Check comments for semantic clues (e.g., "box length", "grid points", "MPI tasks")
3. Examine variable groupings (parameters often appear together)
4. Consider units and value ranges
5. Look for patterns (e.g., X/Y/Z triplets)

RESPONSE FORMAT (JSON):
{{
  "domain_x": ["Lx"],
  "domain_y": ["Ly"],
  "domain_z": ["Lz"],
  "grid_x": ["nxc", "nx"],
  "grid_y": ["nyc"],
  "grid_z": ["nzc"],
  "decomp_x": ["XLEN"],
  "decomp_y": ["YLEN"],
  "decomp_z": ["ZLEN"],
  "timestep": ["dt"],
  "num_cycles": ["ncycles"]
}}

IMPORTANT:
- List ALL variable names that match each canonical type
- If a parameter type is not found, use empty list: []
- Only include variables that actually exist in the file
- Look at BOTH variable names AND comments for meaning
- Consider context (e.g., "number of cells" → grid, not domain)

Respond with ONLY the JSON, no additional text."""
        
        return prompt
    
    def _call_llm(self, prompt: str) -> str:
        """
        Call LLM API with prompt.
        
        Args:
            prompt: Prompt text
            
        Returns:
            LLM response text
        """
        if not HAS_OPENAI:
            raise RuntimeError("OpenAI not available")
        
        try:
            # Try new OpenAI client API (v1.0+)
            from openai import OpenAI  # type: ignore
            client = OpenAI(api_key=self.api_key) if self.api_key else OpenAI()  # type: ignore
            
            response = client.chat.completions.create(  # type: ignore
                model=self.model,
                messages=[
                    {"role": "system", "content": "You are an expert in HPC simulation input files."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.1,  # Low temperature for consistent results
                max_tokens=1000
            )
            
            return response.choices[0].message.content  # type: ignore
            
        except ImportError:
            # Fall back to old API (pre-v1.0)
            import openai as openai_legacy  # type: ignore
            
            if self.api_key:
                openai_legacy.api_key = self.api_key
            
            response = openai_legacy.ChatCompletion.create(  # type: ignore
                model=self.model,
                messages=[
                    {"role": "system", "content": "You are an expert in HPC simulation input files."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.1,
                max_tokens=1000
            )
            
            return response.choices[0].message.content  # type: ignore
    
    def _parse_llm_response(self, response: str) -> Dict[str, List[str]]:
        """
        Parse LLM JSON response into parameter mapping.
        
        Args:
            response: LLM response text
            
        Returns:
            Parameter mapping dictionary
        """
        # Extract JSON from response (handle markdown code blocks)
        json_text = response.strip()
        if json_text.startswith('```'):
            # Remove markdown code fence
            json_text = re.sub(r'```json\s*|\s*```', '', json_text).strip()
        
        try:
            mapping = json.loads(json_text)
            
            # Validate structure
            if not isinstance(mapping, dict):
                raise ValueError("Response is not a dictionary")
            
            # Clean up mapping (remove empty lists, ensure all values are lists)
            cleaned = {}
            for key, value in mapping.items():
                if isinstance(value, list) and len(value) > 0:
                    cleaned[key] = value
                elif isinstance(value, str):
                    cleaned[key] = [value]
            
            return cleaned
            
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse LLM JSON response: {e}")
            logger.debug(f"Response was: {response}")
            raise
    
    def augment_mapping(
        self,
        input_file: Path,
        base_mapping: Dict[str, List[str]]
    ) -> Dict[str, List[str]]:
        """
        Augment existing mapping with LLM-discovered parameters.
        
        This combines rule-based mapping (fast, reliable) with LLM discovery
        (flexible, adaptive) for best results.
        
        Args:
            input_file: Input file to analyze
            base_mapping: Existing parameter mapping
            
        Returns:
            Augmented mapping with LLM-discovered parameters added
        """
        if not self.enabled:
            return base_mapping
        
        logger.info("Augmenting parameter mapping with LLM discovery...")
        
        # Discover new parameters
        llm_mapping = self.discover_parameters(input_file)
        
        # Merge mappings (base takes precedence, LLM adds new ones)
        augmented = base_mapping.copy()
        
        for canonical, variants in llm_mapping.items():
            if canonical in augmented:
                # Add new variants not in base mapping
                existing = set(augmented[canonical])
                new_variants = [v for v in variants if v not in existing]
                if new_variants:
                    augmented[canonical] = augmented[canonical] + new_variants
                    logger.debug(f"Added {len(new_variants)} new variants for {canonical}")
            else:
                # Completely new parameter type
                augmented[canonical] = variants
                logger.info(f"Discovered new parameter type: {canonical}")
        
        return augmented


def create_llm_enhanced_analyzer(
    input_file: Path,
    use_llm: bool = False,
    api_key: Optional[str] = None,
    model: str = "gpt-4",
    yaml_aliases: Optional[Dict[str, List[str]]] = None
):
    """
    Create InputFileAnalyzer with optional LLM enhancement.
    
    Args:
        input_file: Path to input file
        use_llm: Enable LLM-based parameter discovery
        api_key: OpenAI API key (optional)
        model: LLM model to use (default: gpt-4)
        yaml_aliases: Custom parameter mappings from YAML
        
    Returns:
        InputFileAnalyzer instance
    """
    from utils.input_analyzer import InputFileAnalyzer, DEFAULT_PARAMETER_MAP
    
    # Start with default or YAML-provided mapping
    base_mapping = yaml_aliases if yaml_aliases else DEFAULT_PARAMETER_MAP.copy()
    
    if not use_llm:
        # Standard analyzer without LLM (but with YAML aliases if provided)
        return InputFileAnalyzer(input_file, parameter_map=base_mapping)
    
    # LLM-enhanced analyzer
    mapper = LLMParameterMapper(api_key=api_key, model=model)
    
    if mapper.enabled:
        # Use hybrid approach: YAML + LLM
        enhanced_mapping = mapper.discover_with_yaml_hints(
            input_file=input_file,
            yaml_aliases=yaml_aliases
        )
        return InputFileAnalyzer(input_file, parameter_map=enhanced_mapping)
    else:
        # LLM not available, fall back to YAML or standard
        logger.warning("LLM not available - using YAML aliases or default mapping")
        return InputFileAnalyzer(input_file, parameter_map=base_mapping)
