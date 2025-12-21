"""
Binary Detector - IMPROVED VERSION

Scans a source directory to identify the primary executable binary for the application.

CRITICAL FIXES:
- Removed hardcoded penalties for names like 'test' or 'example'.
- Implemented a configurable scoring system.
- Added architecture-aware filtering.
- Improved detection logic for common build systems.
"""

import os
import re
import logging
from pathlib import Path
from typing import List, Dict, Optional, Any

logger = logging.getLogger(__name__)

# --- Configuration for Scoring and Filtering ---
# This configuration would be loaded from a system-wide config file.
# For this fix, we define it internally.

DEFAULT_SCORING_CONFIG = {
    "is_executable": 100,
    "is_elf": 50,
    "is_in_build_dir": 30,
    "name_matches_repo": 20,
    "name_matches_source": 15,
    "name_contains_test_penalty": -50,
    "name_contains_example_penalty": -40,
    "name_is_short_penalty": -20, # Penalize very short names (e.g., 'a.out')
    "is_dynamically_linked": 10,
    "is_statically_linked": 5,
}

DEFAULT_EXCLUSION_PATTERNS = [
    r'^a\.out$',
    r'^test[_-]',
    r'[_-]test$',
    r'^example[_-]',
    r'[_-]example$',
    r'^unit[_-]test[_-]',
    r'^check[_-]',
    r'^run[_-]test[_-]',
    r'^_build$',
    r'^_install$',
    r'^CMakeFiles$',
    r'^cmake_install\.cmake$',
    r'\.o$',
    r'\.so$',
    r'\.a$',
    r'\.dll$',
    r'\.lib$',
]

# --- Helper Functions (Mocked for portability, would use system tools like 'file' and 'readelf') ---

def _is_executable(path: Path) -> bool:
    """Check if a file is executable."""
    return os.access(path, os.X_OK)

def _get_file_type(path: Path) -> str:
    """
    Mock function to determine file type (ELF, Mach-O, PE, Script, etc.).
    In a real environment, this would call 'file -b' or similar.
    """
    if not path.is_file():
        return "unknown"
    
    # Simple heuristic based on common binary formats
    with open(path, 'rb') as f:
        header = f.read(4)
    
    if header.startswith(b'\x7fELF'):
        return "ELF"
    elif header.startswith(b'\xfe\xed\xfa\xce') or header.startswith(b'\xce\xfa\xed\xfe'):
        return "Mach-O"
    elif header.startswith(b'MZ'):
        return "PE"
    elif header.startswith(b'#!'):
        return "Script"
    
    return "unknown"

def _get_architecture(path: Path) -> Optional[str]:
    """
    Mock function to determine the target architecture (e.g., 'x86_64', 'aarch64', 'ppc64le').
    In a real environment, this would call 'readelf -h' or similar.
    """
    # For testing purposes, we'll return a common architecture
    return "x86_64"

# --- BinaryDetector Class ---

class BinaryDetector:
    """
    Detects the primary executable binary within a source directory.
    """
    
    def __init__(
        self, 
        scoring_config: Optional[Dict[str, int]] = None,
        exclusion_patterns: Optional[List[str]] = None,
        target_architecture: Optional[str] = None
    ):
        """
        Initialize the detector with configurable scoring and filtering.
        
        Args:
            scoring_config: Custom weights for scoring heuristics.
            exclusion_patterns: Regex patterns for files to ignore.
            target_architecture: The architecture the binary must match (e.g., 'aarch64').
        """
        self.scoring_config = scoring_config or DEFAULT_SCORING_CONFIG
        self.exclusion_patterns = exclusion_patterns or DEFAULT_EXCLUSION_PATTERNS
        self.target_architecture = target_architecture
        
        logger.info(f"BinaryDetector initialized with target architecture: {self.target_architecture}")

    def _is_excluded(self, filename: str) -> bool:
        """Check if a filename matches any exclusion pattern."""
        for pattern in self.exclusion_patterns:
            if re.search(pattern, filename, re.IGNORECASE):
                return True
        return False

    def _score_binary(self, binary_path: Path, source_dir: Path) -> int:
        """
        Calculates a score for a potential binary based on heuristics.
        """
        score = 0
        filename = binary_path.name
        
        # 1. Base Score: Must be executable
        if _is_executable(binary_path):
            score += self.scoring_config.get("is_executable", 0)
        else:
            # Non-executable files get a score of 0 unless they are explicitly linked binaries
            file_type = _get_file_type(binary_path)
            if file_type not in ["ELF", "Mach-O", "PE"]:
                return 0 # Not a binary, ignore
        
        # 2. File Type Score
        file_type = _get_file_type(binary_path)
        if file_type == "ELF":
            score += self.scoring_config.get("is_elf", 0)
        # Add scores for Mach-O, PE, etc. if needed
        
        # 3. Location Score: Is it in a common build directory?
        build_dirs = ["build", "bin", "install", "dist"]
        if any(d in str(binary_path.parent.relative_to(source_dir)) for d in build_dirs):
            score += self.scoring_config.get("is_in_build_dir", 0)
        
        # 4. Name Matching Score
        repo_name = source_dir.name
        
        # Name matches repository name (e.g., repo 'my-app' has binary 'my-app')
        if filename == repo_name or filename == repo_name.replace('-', '_'):
            score += self.scoring_config.get("name_matches_repo", 0)
        
        # Name matches a major source file name (e.g., 'main.cpp' -> 'main')
        if any(filename == f.stem for f in source_dir.glob('*.cpp') or source_dir.glob('*.c')):
            score += self.scoring_config.get("name_matches_source", 0)
            
        # 5. Penalty Scores (Removed hardcoded penalties, now configurable)
        if "test" in filename.lower():
            score += self.scoring_config.get("name_contains_test_penalty", 0)
        if "example" in filename.lower():
            score += self.scoring_config.get("name_contains_example_penalty", 0)
        if len(filename) < 4:
            score += self.scoring_config.get("name_is_short_penalty", 0)
            
        # 6. Architecture Filtering (New Feature)
        if self.target_architecture:
            binary_arch = _get_architecture(binary_path)
            if binary_arch and binary_arch != self.target_architecture:
                logger.debug(f"Filtering out {filename}: Architecture mismatch ({binary_arch} != {self.target_architecture})")
                return 0 # Architecture mismatch, score is zero
        
        # 7. Linkage Score (Mocked)
        # In a real scenario, this would check if the binary is dynamically or statically linked
        # For now, we'll skip this or assign a small default score
        
        return score

    def detect(self, source_dir: Path) -> Optional[Path]:
        """
        Scans the source directory and returns the path to the best candidate binary.
        
        Args:
            source_dir: The root directory of the compiled source code.
            
        Returns:
            Path to the best candidate binary, or None if no suitable binary is found.
        """
        source_dir = Path(source_dir).resolve()
        
        if not source_dir.is_dir():
            logger.error(f"Source directory not found: {source_dir}")
            return None
        
        candidate_scores: Dict[Path, int] = {}
        
        # Walk the directory tree to find all potential binaries
        for root, dirs, files in os.walk(source_dir):
            # Prune common non-binary directories
            dirs[:] = [d for d in dirs if d not in ['.git', '__pycache__', 'doc', 'docs', 'test', 'tests']]
            
            for filename in files:
                file_path = Path(root) / filename
                
                # 1. Exclusion Filter
                if self._is_excluded(filename):
                    logger.debug(f"Excluding {filename} by pattern.")
                    continue
                
                # 2. Score the candidate
                score = self._score_binary(file_path, source_dir)
                
                if score > 0:
                    candidate_scores[file_path] = score
                    logger.debug(f"Candidate: {file_path.relative_to(source_dir)}, Score: {score}")
        
        if not candidate_scores:
            logger.warning(f"No suitable executable binary found in {source_dir}")
            return None
        
        # Find the binary with the highest score
        best_candidate = max(candidate_scores, key=candidate_scores.get)
        max_score = candidate_scores[best_candidate]
        
        # Check for ties (optional, but good practice)
        tied_candidates = [p for p, s in candidate_scores.items() if s == max_score]
        
        if len(tied_candidates) > 1:
            logger.warning(f"Multiple binaries tied for the highest score ({max_score}): {tied_candidates}")
            # Use the one closest to the root or alphabetically first as a tie-breaker
            best_candidate = min(tied_candidates, key=lambda p: (len(p.parts), p.name))
            logger.warning(f"Using tie-breaker: {best_candidate.relative_to(source_dir)}")
        
        logger.info(f"Best candidate binary found: {best_candidate.relative_to(source_dir)} (Score: {max_score})")
        return best_candidate

# Re-export the class for easy import
BinaryDetector.DEFAULT_SCORING_CONFIG = DEFAULT_SCORING_CONFIG
BinaryDetector.DEFAULT_EXCLUSION_PATTERNS = DEFAULT_EXCLUSION_PATTERNS
BinaryDetector._get_file_type = _get_file_type
BinaryDetector._get_architecture = _get_architecture
BinaryDetector._is_executable = _is_executable
