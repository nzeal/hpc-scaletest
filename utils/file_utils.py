"""
File utils.
"""

import shutil
import logging
from pathlib import Path
from typing import Optional


logger = logging.getLogger(__name__)


def create_directory(path: Path, exist_ok: bool = True) -> bool:
    """
    Create dir.
    """
    try:
        path.mkdir(parents=True, exist_ok=exist_ok)
        logger.debug(f"Created dir: {path}")
        return True
    except Exception as e:
        logger.error(f"Failed to create dir {path}: {e}")
        return False


def write_file(path: Path, content: str) -> bool:
    """
    Write file.
    """
    try:
        path.write_text(content)
        logger.debug(f"Wrote file: {path}")
        return True
    except Exception as e:
        logger.error(f"Failed to write file {path}: {e}")
        return False
