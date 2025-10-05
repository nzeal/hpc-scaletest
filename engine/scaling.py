# engine/scaling.py
import math
from typing import Dict, Tuple, Any
from core.types import ScalingType
from core.config import ScalingConfig

class ScalingEngine:
    """
    Engine for generating scaling configurations for strong and weak scaling tests.
    Supports 3D domain decomposition with alternating dimension doubling.
    """

    def generate_scaling_configs(
        self, 
        config: ScalingConfig, 
        procs_per_node: int = 128
    ) -> Dict[int, Dict[str, Any]]:
        """
        Generate scaling configurations for node counts up to max_nodes.

        :param config: ScalingConfig with initial parameters
        :param procs_per_node: Processes per node for node calculation
        :return: Dict of {actual_nodes: config_dict}
        """
        configs = {}
        # Initial values
        xproc, yproc, zproc = config.initial_procs
        lx, ly, lz = config.initial_domain
        nxc, nyc, nzc = config.initial_cells
        
        # Generate node list: 1, 2, 4, ..., up to <= max_nodes
        max_log = int(math.log2(config.max_nodes))
        nodes_list = [1] + [2**i for i in range(1, max_log + 1) if 2**i <= config.max_nodes]

        initial_total_procs = xproc * yproc * zproc
        initial_nodes = math.ceil(initial_total_procs / procs_per_node)
        if initial_nodes > 1:
            nodes_list = [initial_nodes] + nodes_list  # Adjust if initial >1
        nodes_list = sorted(set(nodes_list))  # Unique sorted

        current_x, current_y, current_z = xproc, yproc, zproc
        current_lx, current_ly, current_lz = lx, ly, lz
        current_nxc, current_nyc, current_nzc = nxc, nyc, nzc

        for j, target_nodes in enumerate(nodes_list):
            # Calculate required total_procs for target_nodes
            target_total_procs = target_nodes * procs_per_node
            # Scale procs from initial, preserving aspect ratio roughly
            scale_procs = target_total_procs / initial_total_procs
            # Alternate doubling: for simplicity, double x then y repeatedly
            while current_x * current_y * current_z < target_total_procs:
                if (j % 2 == 0):
                    current_x = min(current_x * 2, target_total_procs // (current_y * current_z))
                else:
                    current_y = min(current_y * 2, target_total_procs // (current_x * current_z))
                # Ensure exact
                current_x = max(1, target_total_procs // (current_y * current_z))
            
            actual_total_procs = current_x * current_y * current_z
            actual_nodes = math.ceil(actual_total_procs / procs_per_node)

            # Adjust sizes based on scaling type
            if config.scaling_type == ScalingType.STRONG:
                # Constant problem size
                pass
            elif config.scaling_type == ScalingType.WEAK:
                # Scale problem size proportionally to procs
                scale_factor = (actual_total_procs / initial_total_procs) ** (1/3)  # Cubic root for 3D
                current_lx = lx * scale_factor
                current_ly = ly * scale_factor
                current_lz = lz * scale_factor
                current_nxc = int(nxc * scale_factor)
                current_nyc = int(nyc * scale_factor)
                current_nzc = int(nzc * scale_factor)

            configs[actual_nodes] = {
                'n_nodes': actual_nodes,
                'xproc': current_x,
                'yproc': current_y,
                'zproc': current_z,
                'lx': current_lx,
                'ly': current_ly,
                'lz': current_lz,
                'nxc': current_nxc,
                'nyc': current_nyc,
                'nzc': current_nzc,
                'total_procs': actual_total_procs
            }

            # Reset for next if needed, but since progressive, carry over

        return configs
