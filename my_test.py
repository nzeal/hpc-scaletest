from pathlib import Path
from core.test_definition import Test

def create_test():
    test = Test(
        name="my_benchmark",
        input_file=Path("./input.dat"),  # Your benchmark input file
        command=["./my_app"]  # Executable and base args
    )
    
    test.set_backend(
        scheduler="slurm",  # Or "local"
        launcher="srun",    # Or "mpirun"
        module_system="lmod"
    )
    
    test.set_resources(
        max_nodes=64,
        procs_per_node=128,
        time_limit="01:00:00",
        partition="gpu"  # If using GPUs
    )
    
    test.set_scaling(
        scaling_type="strong",  # Or "weak"
        max_nodes=64,
        initial_procs=(2, 2, 2),
        initial_domain=(10.0, 10.0, 10.0),
        initial_cells=(256, 256, 256)
    )
    
    test.set_modules(["gcc/11.2.0", "openmpi/4.1.1"])
    test.set_env({"OMP_NUM_THREADS": "1"})
    
    return test
