"""
TCL-based Environment Modules 4.x backend.
"""

import subprocess
from typing import List, Optional

from core.abstracts import ModuleSystemInterface


class TMod4Backend(ModuleSystemInterface):
    def generate_load_commands(self, modules: List[str]) -> List[str]:
        return [f"module load {module}" for module in modules]
    
    def generate_unload_commands(self, modules: List[str]) -> List[str]:
        return [f"module unload {module}" for module in modules]
    
    def list_available_modules(self, pattern: Optional[str] = None) -> List[str]:
        try:
            cmd = ["module", "-t", "avail"]
            if pattern:
                cmd.append(pattern)
            result = subprocess.run(cmd, capture_output=True, text=True, shell=False)
            lines = result.stderr.split('\n')
            modules = []
            for line in lines:
                line = line.strip()
                if line and not line.startswith('-') and not line.endswith(':'):
                    modules.append(line)
            return modules
        except Exception:
            return []
    
    def is_module_available(self, module: str) -> bool:
        try:
            result = subprocess.run(["module", "is-avail", module], capture_output=True)
            return result.returncode == 0
        except Exception:
            return False
