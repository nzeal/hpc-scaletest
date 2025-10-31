"""
No-op module system backend.
"""

from typing import List, Optional

from core.abstracts import ModuleSystemInterface


class NoModBackend(ModuleSystemInterface):
    def generate_load_commands(self, modules: List[str]) -> List[str]:
        return []
    
    def generate_unload_commands(self, modules: List[str]) -> List[str]:
        return []
    
    def list_available_modules(self, pattern: Optional[str] = None) -> List[str]:
        return []
    
    def is_module_available(self, module: str) -> bool:
        return False
