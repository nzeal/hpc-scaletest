# backends/modules/nomod.py
from typing import List
from ...core.abstracts import AbstractModuleSystem
import logging

logger = logging.getLogger(__name__)

class NoModBackend(AbstractModuleSystem):
    """No-op backend for systems without module support."""

    def load(self, modules: List[str]) -> None:
        """Log load action without executing."""
        if modules:
            logger.info(f"NoMod: Would load modules: {', '.join(modules)}")
        else:
            logger.info("NoMod: No modules to load")

    def unload(self, modules: List[str]) -> None:
        """Log unload action without executing."""
        if modules:
            logger.info(f"NoMod: Would unload modules: {', '.join(modules)}")
        else:
            logger.info("NoMod: No modules to unload")

    def purge(self) -> None:
        """Log purge action without executing."""
        logger.info("NoMod: Would purge all modules")

    def list(self) -> List[str]:
        """Return empty list as no modules are managed."""
        logger.debug("NoMod: Listing modules (empty)")
        return []
