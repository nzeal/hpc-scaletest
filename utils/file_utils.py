# utils/file_utils.py
from pathlib import Path
from typing import List, Tuple

def modify_file(file_path: Path, replacements: List[Tuple[str, str]]) -> None:
    """
    Modify a file by replacing specified lines or patterns (sed-like).

    :param file_path: Path to the file
    :param replacements: List of (old_text, new_text) tuples
    """
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    content = file_path.read_text(encoding='utf-8')
    for old, new in replacements:
        content = content.replace(old, new)
    
    file_path.write_text(content, encoding='utf-8')
    logger.debug(f"Modified file {file_path} with {len(replacements)} replacements")
