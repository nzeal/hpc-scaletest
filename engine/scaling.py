"""
Scaling generator.
"""

import logging
from typing import List

from core.config import ScalingConfig, JobConfig, ResourceConfig
from core.types import ScalingType, ProcsDecomposition, DomainSize, CellCount


logger = logging.getLogger(__name__)


class ScalingEngine:
    def __init__(self, scaling_config: ScalingConfig, resource_config: ResourceConfig):
        self.scaling_config = scaling_config
        self.resource_config = resource_config
    
    def generate_job_configs(self) -> List[JobConfig]:
        node_sequence = self.scaling_config.get_node_sequence()
        if self.scaling_config.scaling_type == ScalingType.STRONG:
            return self._generate_strong_scaling(node_sequence)
        return self._generate_weak_scaling(node_sequence)
    
    def _generate_strong_scaling(self, node_sequence: List[int]) -> List[JobConfig]:
        configs = []
        base_px, base_py, base_pz = self.scaling_config.initial_procs
        for num_nodes in node_sequence:
            # Start from initial decomposition for each node count
            px, py, pz = base_px, base_py, base_pz
            total_procs = num_nodes * self.resource_config.procs_per_node
            current_procs = px * py * pz
            while current_procs < total_procs:
                if px <= py:
                    px *= 2
                else:
                    py *= 2
                current_procs = px * py * pz
            configs.append(JobConfig(
                job_id=f"nodes_{num_nodes}",
                num_nodes=num_nodes,
                num_procs=total_procs,
                procs_decomposition=(px, py, pz),
                domain_size=self.scaling_config.initial_domain,
                cell_count=self.scaling_config.initial_cells
            ))
            logger.debug(f"Strong: {num_nodes} nodes, procs=({px},{py},{pz})")
        return configs
    
    def _generate_weak_scaling(self, node_sequence: List[int]) -> List[JobConfig]:
        configs = []
        base_px, base_py, base_pz = self.scaling_config.initial_procs
        initial_domain = self.scaling_config.initial_domain
        initial_cells = self.scaling_config.initial_cells
        initial_total = base_px * base_py * base_pz
        for num_nodes in node_sequence:
            # Start from initial decomposition for each node count
            px, py, pz = base_px, base_py, base_pz
            total_procs = num_nodes * self.resource_config.procs_per_node
            scale_factor = total_procs / initial_total
            current_procs = px * py * pz
            while current_procs < total_procs:
                if px <= py:
                    px *= 2
                else:
                    py *= 2
                current_procs = px * py * pz
            scaled_domain = None
            scaled_cells = None
            if initial_domain:
                dx, dy, dz = initial_domain
                scaled_domain = (dx * scale_factor, dy * scale_factor, dz * scale_factor)
            if initial_cells:
                nx, ny, nz = initial_cells
                scaled_cells = (int(nx * scale_factor), int(ny * scale_factor), int(nz * scale_factor))
            configs.append(JobConfig(
                job_id=f"nodes_{num_nodes}",
                num_nodes=num_nodes,
                num_procs=total_procs,
                procs_decomposition=(px, py, pz),
                domain_size=scaled_domain,
                cell_count=scaled_cells
            ))
            logger.debug(f"Weak: {num_nodes} nodes, procs=({px},{py},{pz})")
        return configs
