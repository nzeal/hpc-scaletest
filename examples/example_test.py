"""
Example test definition for HPC-ScaleTest.

This demonstrates how users define benchmark tests without worrying
about scheduler details or backend implementations.
"""

# examples/example_test.py
from pathlib import Path
from core.test_definition import Test

def create_test():
    test = Test(
        name="example_benchmark",
        input_file=Path("./input.dat"),
        command=["./example_app"]
    )
    
    test.set_backend(
        scheduler="local",  # For example, use local for testing
        launcher="mpirun",
        module_system="nomod"
    )
    
    test.set_resources(
        max_nodes=4,  # Small for example
        procs_per_node=4,
        gpus_per_node=0,
        time_limit="00:10:00"
    )
    
    test.set_scaling(
        scaling_type="strong",
        max_nodes=4,
        initial_procs=(1, 1, 1),
        initial_domain=(10.0, 10.0, 10.0),
        initial_cells=(256, 256, 256)
    )
    
    test.set_modules(["gcc/11.2.0"])  # Example modules
    test.set_env({"OMP_NUM_THREADS": "1"})
    
    return test

if __name__ == "__main__":
    test = create_test()
    from engine.runner import TestRunner
    runner = TestRunner(test)
    runner.run()
