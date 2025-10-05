from pathlib import Path
from core.test_definition import Test

# Create test
test = Test(
    name="hello_mpi",
    input_file=Path("input.txt"),
    command=["mpirun", "-n", "4", "hostname"]
)

# Configure for local testing
test.set_backend(scheduler="local", launcher="mpirun")
test.set_resources(max_nodes=1, procs_per_node=4)
test.set_scaling(
    scaling_type="strong",
    max_nodes=1,
    initial_procs=(1, 1, 1)
)
