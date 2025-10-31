#!/usr/bin/env python3
"""
Advanced workflow example for HPC Auto.
Demonstrates full control over the orchestrator.
"""

from pathlib import Path
from engine.orchestrator import HPCOrchestrator, OrchestratorConfig
from utils.logging_config import setup_logging
import logging


def main():
    """Run an advanced scaling test with full configuration."""
    # Setup logging
    setup_logging(level=logging.INFO)
    logger = logging.getLogger(__name__)
    
    # Create detailed configuration
    config = OrchestratorConfig(
        # Source configuration
        source="/path/to/your/code",
        
        # Test configuration
        scaling_type="strong",
        hardware_type="cpu",
        test_name="my_advanced_test",
        
        # Scaling parameters
        max_nodes=32,
        initial_procs=(4, 4, 2),  # 3D decomposition
        initial_domain=(40.96, 20.48, 1.0),  # Domain size
        initial_cells=(896, 512, 1),  # Cell count
        
        # Resource configuration
        procs_per_node=128,
        gpus_per_node=0,
        time_limit="04:00:00",
        partition="production",
        account="my_project",
        
        # Backend configuration
        scheduler="slurm",
        launcher="srun",
        module_system="lmod",
        
        # Build configuration
        build_system="cmake",
        modules=[
            "gcc/11.2.0",
            "openmpi/4.1.1",
            "hdf5/1.12.0"
        ],
        build_flags={
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_C_COMPILER": "mpicc",
            "CMAKE_CXX_COMPILER": "mpicxx",
            "ENABLE_HDF5": "ON"
        },
        
        # Behavior configuration
        auto_submit_jobs=False,  # Generate scripts but don't submit
        cleanup_after_build=False,
        generate_reports=True,
        
        # Output configuration
        output_dir=Path("results/advanced_test"),
        workspace_dir=Path("workspace")
    )
    
    logger.info("Configuration created")
    logger.info(f"  Test: {config.test_name}")
    logger.info(f"  Scaling: {config.scaling_type}")
    logger.info(f"  Max nodes: {config.max_nodes}")
    logger.info(f"  Auto-submit: {config.auto_submit_jobs}")
    
    # Create orchestrator
    orchestrator = HPCOrchestrator(config)
    
    # Run workflow
    success = orchestrator.run()
    
    if success:
        logger.info("\n" + "="*60)
        logger.info("WORKFLOW COMPLETED SUCCESSFULLY")
        logger.info("="*60)
        logger.info(f"Results saved to: {config.output_dir}")
        
        if not config.auto_submit_jobs:
            logger.info("\nJob scripts generated but not submitted.")
            logger.info("To submit jobs manually:")
            logger.info(f"  python utils/job_submitter.py {config.output_dir}")
    else:
        logger.error("Workflow failed!")
        return 1
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
