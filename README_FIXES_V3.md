# HPC-ScaleTest Fixes - Version 3

## Summary of Issues Fixed

### 1. Launcher Not Respected (srun vs mpirun)
**Problem**: YAML `launcher: srun` was ignored, always using mpirun.
**Solution**: Fixed code flow to pass launcher setting to scheduler.

### 2. CPU Jobs Using Wrong MPI Options  
**Problem**: CPU jobs got `--map-by ppr:112:node:PE=1` which shouldn't be used.
**Solution**: Only GPU jobs with mpirun get `--map-by` options.

### 3. GPU bind.sh Not Created
**Problem**: bind.sh was not created, causing GPU jobs to fail.
**Solution**: GPU wrapper script (`gpu_wrapper.sh`) is created inline in job.sh.

### 4. QOSMinCpuNotSatisfied Error
**Problem**: SLURM rejected jobs due to insufficient CPU allocation.
**Solution**: Changed SLURM directives to allocate full node CPUs.

## Generated Commands

### CPU Jobs
```bash
# With srun (recommended for SLURM):
srun --ntasks=112 --ntasks-per-node=112 $BINARY/iPIC3D os-stdin

# With mpirun:
mpirun -np 112 $BINARY/iPIC3D os-stdin
```

### GPU Jobs  
```bash
# With mpirun (includes gpu_wrapper.sh for GPU binding):
mpirun -np 4 --map-by ppr:4:node:PE=8 --report-bindings ./gpu_wrapper.sh $BINARY/iPIC3D os-stdin

# With srun:
srun --ntasks=4 --ntasks-per-node=4 --cpus-per-task=8 $BINARY/iPIC3D os-stdin
```

## SLURM Directives Generated

### CPU Job (112 cores):
```bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=112
#SBATCH --cpus-per-task=1
```

### GPU Job (32 cores, 4 GPUs):
```bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4
```

## YAML Configuration Examples

### CPU Weak Scaling with srun:
```yaml
repository: https://github.com/user/iPIC3D-CPU.git
scaling: weak
nodes: 32
launcher: srun           # <-- Use SLURM's native launcher
hardware: cpu
procs_per_node: 112
partition: dcgp_usr_prod
account: cin_staff
```

### GPU Weak Scaling with mpirun:
```yaml
repository: https://github.com/user/iPIC3D-GPU.git
scaling: weak
nodes: 16
launcher: mpirun        # <-- Use OpenMPI
hardware: gpu
gpus_per_node: 4
procs_per_node: 32      # <-- IMPORTANT: Set to total CPU cores per node
partition: boost_usr_prod
account: cin_staff
```

## GPU Binding (Inline in job.sh)

For GPU jobs with mpirun, the job script creates `gpu_wrapper.sh` at runtime:
```bash
cat > gpu_wrapper.sh << 'WRAPPER_EOF'
#!/bin/bash
if [ -n "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
    LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
elif [ -n "$MPI_LOCALRANKID" ]; then
    LOCAL_RANK=$MPI_LOCALRANKID
else
    LOCAL_RANK=0
fi
export CUDA_VISIBLE_DEVICES=$LOCAL_RANK
export ROCR_VISIBLE_DEVICES=$LOCAL_RANK
export HIP_VISIBLE_DEVICES=$LOCAL_RANK
exec "$@"
WRAPPER_EOF
chmod +x gpu_wrapper.sh
```

## Testing

After extracting the zip, test with:
```bash
cd hpc-scaletest-fixed
python3 -c "
from backends.schedulers.slurm import SlurmScheduler
from core.config import JobConfig, ResourceConfig
from pathlib import Path

# Test CPU with srun
scheduler = SlurmScheduler({'launcher': 'srun'})
resource_config = ResourceConfig()
resource_config.procs_per_node = 112
resource_config.gpus_per_node = 0
job_config = JobConfig(job_id='test', num_nodes=1, num_procs=112, procs_decomposition=(8,7,2), working_dir=Path('/tmp'))
script = scheduler.generate_job_script(job_config, resource_config, ['\$BINARY/iPIC3D', 'os-stdin'], [])
print('Generated script (showing run command):')
for line in script.split('\n'):
    if 'srun' in line and 'echo' not in line:
        print(line)
"
```

Expected output:
```
# Launcher: srun
srun --ntasks=112 --ntasks-per-node=112 $BINARY/iPIC3D os-stdin
```

## Files Modified

1. **backends/schedulers/slurm.py**
   - Added launcher option to constructor
   - Rewrote generate_job_script() to:
     - Respect launcher setting
     - Generate correct commands per job type
     - Create GPU wrapper inline
     - Use proper SLURM directives

2. **engine/runner.py**
   - Pass launcher to scheduler options

## Troubleshooting

### Still getting QOSMinCpuNotSatisfied?
Check that `procs_per_node` in your YAML matches your partition's CPU count:
- Leonardo DCGP: `procs_per_node: 112`
- Leonardo Booster: `procs_per_node: 32`

### Jobs failing with "bind.sh not found"?
The latest code creates `gpu_wrapper.sh` inline. If you see this error, you're running old code.
