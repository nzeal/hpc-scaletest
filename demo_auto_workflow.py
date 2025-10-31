#!/usr/bin/env python3
"""
Demo script for the automated HPC scaling workflow.
Shows how the framework can be used with minimal configuration.
"""

import sys
from pathlib import Path
from engine.orchestrator import create_simple_workflow, HPCOrchestrator, OrchestratorConfig
from utils.logging_config import setup_logging
import logging


def demo_simple_workflow():
    """Demonstrate the simplest possible workflow."""
    print("\n" + "="*70)
    print("DEMO 1: Simple Workflow - Minimal Configuration")
    print("="*70)
    
    print("\nThis demo shows how to run a scaling test with minimal input.")
    print("The framework will:")
    print("  1. Analyze the code and detect dependencies")
    print("  2. Configure the build system")
    print("  3. Generate scaling test configurations")
    print("  4. Create job submission scripts")
    print("  5. Generate efficiency reports")
    
    print("\nExample usage:")
    print("  orchestrator = create_simple_workflow(")
    print("      source='/path/to/code',")
    print("      scaling_type='strong',")
    print("      max_nodes=8")
    print("  )")
    print("  orchestrator.run()")
    
    print("\nThis would:")
    print("  - Auto-detect build system from project files")
    print("  - Auto-detect dependencies from README")
    print("  - Generate jobs for 1, 2, 4, 8 nodes")
    print("  - Submit jobs and monitor completion")
    print("  - Generate efficiency report")


def demo_advanced_workflow():
    """Demonstrate advanced workflow with full configuration."""
    print("\n" + "="*70)
    print("DEMO 2: Advanced Workflow - Full Configuration")
    print("="*70)
    
    print("\nThis demo shows full control over the workflow.")
    
    config_example = """
config = OrchestratorConfig(
    # Source
    source="/path/to/code",
    
    # Test configuration
    scaling_type="strong",
    hardware_type="cpu",
    test_name="my_benchmark",
    
    # Scaling parameters
    max_nodes=32,
    initial_procs=(4, 4, 2),
    initial_domain=(40.96, 20.48, 1.0),
    initial_cells=(896, 512, 1),
    
    # Resources
    procs_per_node=128,
    time_limit="04:00:00",
    partition="production",
    account="my_project",
    
    # Backend
    scheduler="slurm",
    launcher="srun",
    module_system="lmod",
    
    # Build
    build_system="cmake",
    modules=["gcc/11.2.0", "openmpi/4.1.1"],
    build_flags={
        "CMAKE_BUILD_TYPE": "Release",
        "ENABLE_HDF5": "ON"
    },
    
    # Behavior
    auto_submit_jobs=True,
    generate_reports=True
)

orchestrator = HPCOrchestrator(config)
orchestrator.run()
"""
    
    print("\nExample configuration:")
    print(config_example)


def demo_cli_usage():
    """Demonstrate CLI usage."""
    print("\n" + "="*70)
    print("DEMO 3: Command-Line Interface")
    print("="*70)
    
    print("\nThe CLI provides the simplest interface:")
    
    examples = [
        ("Basic strong scaling", 
         "python hpc_auto.py /path/to/code --scaling strong --nodes 8"),
        
        ("Clone from Git",
         "python hpc_auto.py https://github.com/user/app.git --scaling weak --nodes 16"),
        
        ("GPU test",
         "python hpc_auto.py /path/to/cuda-app --hardware gpu --nodes 4"),
        
        ("Custom resources",
         "python hpc_auto.py /path/to/code --nodes 16 --procs-per-node 64 --partition gpu"),
        
        ("Manual modules",
         "python hpc_auto.py /path/to/code --modules 'gcc/11.2.0,openmpi/4.1.1' --nodes 8"),
        
        ("Generate scripts only",
         "python hpc_auto.py /path/to/code --nodes 4 --no-submit"),
    ]
    
    for title, command in examples:
        print(f"\n{title}:")
        print(f"  $ {command}")


def demo_workflow_steps():
    """Demonstrate the workflow steps."""
    print("\n" + "="*70)
    print("DEMO 4: Workflow Steps Breakdown")
    print("="*70)
    
    steps = [
        ("Step 1: Code Acquisition", [
            "- Accepts local path or Git URL",
            "- Clones Git repos to workspace/",
            "- Validates source availability"
        ]),
        
        ("Step 2: Code Analysis", [
            "- Scans README.md and documentation",
            "- Detects build system (CMake, Make, etc.)",
            "- Identifies MPI, CUDA, OpenMP requirements",
            "- Extracts module dependencies",
            "- Confidence score: 0.0 - 1.0"
        ]),
        
        ("Step 3: Compilation", [
            "- Loads detected modules",
            "- Configures build with proper flags",
            "- Compiles with detected build system",
            "- Installs to output directory"
        ]),
        
        ("Step 4: Test Generation", [
            "- Creates node sequence (1, 2, 4, 8, ...)",
            "- Generates process decompositions",
            "- Scales input files (for weak scaling)",
            "- Creates job submission scripts"
        ]),
        
        ("Step 5: Job Execution", [
            "- Submits jobs to scheduler",
            "- Monitors job status",
            "- Collects output files",
            "- Parses timing data"
        ]),
        
        ("Step 6: Report Generation", [
            "- Calculates speedup metrics",
            "- Computes efficiency percentages",
            "- Generates text report",
            "- Creates JSON data file"
        ])
    ]
    
    for step_title, details in steps:
        print(f"\n{step_title}:")
        for detail in details:
            print(f"  {detail}")


def demo_output_structure():
    """Show the output structure."""
    print("\n" + "="*70)
    print("DEMO 5: Output Structure")
    print("="*70)
    
    structure = """
output/
└── my_app_strong_20251025_104500/
    ├── nodes_1/
    │   ├── job.sh                  # Slurm submission script
    │   ├── input.dat               # Application input (scaled)
    │   ├── metadata.json           # Job configuration
    │   └── out_*.out               # Job output (after execution)
    │
    ├── nodes_2/
    │   └── (same structure)
    │
    ├── nodes_4/
    │   └── (same structure)
    │
    ├── nodes_8/
    │   └── (same structure)
    │
    ├── install/                    # Compiled executable
    │   └── my_app
    │
    ├── summary.json                # Overall test summary
    ├── StrongScalingReport.txt     # Human-readable report
    └── strong_scaling_report.json  # Machine-readable data
"""
    
    print(structure)
    
    print("\nExample StrongScalingReport.txt:")
    report = """
==============================================================================
Strong Scaling Efficiency Report
==============================================================================
Test Name: my_app
Generated: 2025-10-25 10:45:00
==============================================================================

Nodes      Procs      Time(s)      Speedup      Efficiency   Decomposition
------------------------------------------------------------------------------
1          128        100.00       1.00         100.0%       (8×8×2)
2          256        52.00        1.92         96.2%        (16×8×2)
4          512        28.00        3.57         89.3%        (16×16×2)
8          1024       16.00        6.25         78.1%        (32×16×2)
==============================================================================

Summary Statistics:
----------------------------------------
  Average Efficiency: 90.9%
  Maximum Efficiency: 100.0%
  Minimum Efficiency: 78.1%
  Maximum Speedup: 6.25x
==============================================================================
"""
    print(report)


def demo_comparison():
    """Compare old vs new approach."""
    print("\n" + "="*70)
    print("DEMO 6: Old vs New Approach Comparison")
    print("="*70)
    
    print("\n### OLD APPROACH (Manual) ###")
    print("1. Clone/copy code manually")
    print("2. Read README to find dependencies")
    print("3. Load modules manually")
    print("4. Configure build system manually")
    print("5. Compile code")
    print("6. Write Test definition in Python")
    print("7. Configure backend, resources, scaling")
    print("8. Run test")
    print("9. Wait for completion")
    print("10. Manually parse outputs")
    print("11. Manually calculate metrics")
    print("12. Create report in spreadsheet")
    print("\n⏱️  Time: 2-4 hours")
    print("❌ Error-prone")
    print("❌ Not reproducible")
    
    print("\n### NEW APPROACH (Automated) ###")
    print("$ python hpc_auto.py /path/to/code --scaling strong --nodes 8")
    print("\n✅ All steps automated")
    print("✅ Intelligent dependency detection")
    print("✅ Beautiful efficiency reports")
    print("⏱️  Time: Single command, ~5 minutes setup")
    print("✨ Fully reproducible")


def main():
    """Run all demos."""
    setup_logging(level=logging.INFO)
    
    print("\n" + "="*70)
    print("HPC AUTO - AUTOMATED SCALING FRAMEWORK DEMO")
    print("="*70)
    print("\nThis demo showcases the new automated workflow capabilities.")
    print("No actual tests will be run - this is just a demonstration.")
    
    demos = [
        demo_simple_workflow,
        demo_advanced_workflow,
        demo_cli_usage,
        demo_workflow_steps,
        demo_output_structure,
        demo_comparison
    ]
    
    for demo_func in demos:
        demo_func()
        input("\nPress Enter to continue to next demo...")
    
    print("\n" + "="*70)
    print("DEMO COMPLETE")
    print("="*70)
    print("\nTo get started with real tests:")
    print("  1. Read: docs/QUICKSTART.md")
    print("  2. Read: docs/AUTOMATED_WORKFLOW.md")
    print("  3. Run: python hpc_auto.py --help")
    print("  4. Try: python hpc_auto.py /path/to/code --nodes 2 --no-submit")
    print("\n" + "="*70)


if __name__ == '__main__':
    main()
