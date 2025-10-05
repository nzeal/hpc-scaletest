# tests/test_scaling.py
import unittest
import math
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from engine.scaling import ScalingEngine
from core.types import ScalingConfig, ScalingType

class TestScaling(unittest.TestCase):
    def setUp(self):
        self.engine = ScalingEngine()
        self.base_config = ScalingConfig(
            scaling_type=ScalingType.STRONG,
            max_nodes=8,
            initial_procs=(1, 1, 1),
            initial_domain=(10.0, 10.0, 10.0),
            initial_cells=(256, 256, 256)
        )

    def test_generate_configs_strong(self):
        # Use small procs_per_node for predictable totals
        configs = self.engine.generate_scaling_configs(self.base_config, procs_per_node=4)
        self.assertEqual(len(configs), 4)
        self.assertIn(1, configs)
        self.assertEqual(configs[1]['n_nodes'], 1)
        self.assertEqual(configs[1]['total_procs'], 4)
        self.assertEqual(configs[1]['xproc'], 4)
        self.assertEqual(configs[1]['yproc'], 1)
        self.assertEqual(configs[1]['zproc'], 1)
        self.assertEqual(configs[1]['lx'], 10.0)  # Constant for strong
        self.assertEqual(configs[1]['nxc'], 256)

        # Check progressive doubling
        self.assertEqual(configs[2]['total_procs'], 8)
        self.assertEqual(configs[4]['total_procs'], 16)
        self.assertEqual(configs[8]['total_procs'], 32)

    def test_generate_configs_weak(self):
        weak_config = ScalingConfig(
            scaling_type=ScalingType.WEAK,
            max_nodes=8,
            initial_procs=(1, 1, 1),
            initial_domain=(10.0, 10.0, 10.0),
            initial_cells=(256, 256, 256)
        )
        configs = self.engine.generate_scaling_configs(weak_config, procs_per_node=4)
        self.assertEqual(len(configs), 4)
        # Check scaling: scale_factor = (total_procs / 1)**(1/3)
        # For total=4: ~1.587, lx~15.87, nxc~406
        self.assertAlmostEqual(configs[1]['lx'], 15.874, places=2)
        self.assertEqual(configs[1]['nxc'], 406)
        # For total=8: 2.0, lx=20.0, nxc=512
        self.assertEqual(configs[2]['lx'], 20.0)
        self.assertEqual(configs[2]['nxc'], 512)

    def test_initial_procs_greater_than_one(self):
        config = ScalingConfig(
            scaling_type=ScalingType.STRONG,
            max_nodes=8,
            initial_procs=(2, 2, 1),  # total=4
            initial_domain=(10.0, 10.0, 10.0),
            initial_cells=(256, 256, 256)
        )
        configs = self.engine.generate_scaling_configs(config, procs_per_node=4)
        # initial_nodes=1 (ceil(4/4)=1), but starts with total=4
        self.assertEqual(list(configs.keys())[0], 1)
        self.assertEqual(configs[1]['total_procs'], 4)
        self.assertEqual(configs[1]['xproc'], 2)
        self.assertEqual(configs[1]['yproc'], 2)

    def test_max_nodes_one(self):
        config = ScalingConfig(
            scaling_type=ScalingType.STRONG,
            max_nodes=1,
            initial_procs=(1, 1, 1),
            initial_domain=(10.0, 10.0, 10.0),
            initial_cells=(256, 256, 256)
        )
        configs = self.engine.generate_scaling_configs(config, procs_per_node=4)
        self.assertEqual(len(configs), 1)
        self.assertEqual(configs[1]['total_procs'], 4)

if __name__ == '__main__':
    unittest.main()
