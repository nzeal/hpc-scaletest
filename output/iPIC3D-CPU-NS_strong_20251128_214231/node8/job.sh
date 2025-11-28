#!/bin/bash

#SBATCH --nodes=8
#SBATCH --partition=dcgp_usr_prod
#SBATCH --qos=normal
#SBATCH --ntasks=896
#SBATCH --ntasks-per-node=112
#SBATCH --cpus-per-task=1
#SBATCH --gres=tmpfs:10g
#SBATCH -A cin_staff
#SBATCH --time=02:00:00
#SBATCH --job-name=node8
#SBATCH -o /leonardo_scratch/large/userinternal/nshukla1/finalhpc/hpc-scaletest/output/iPIC3D-CPU-NS_strong_20251128_214231/node8/job.out
#SBATCH -e /leonardo_scratch/large/userinternal/nshukla1/finalhpc/hpc-scaletest/output/iPIC3D-CPU-NS_strong_20251128_214231/node8/job.err

# Environment setup
module load hdf5/1.14.3--intel-oneapi-mpi--2021.12.1--oneapi--2024.1.0
export OMP_NUM_THREADS=1

echo "Loading modules"

# Change to working directory
cd /leonardo_scratch/large/userinternal/nshukla1/finalhpc/hpc-scaletest/output/iPIC3D-CPU-NS_strong_20251128_214231/node8
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Binary location
BINARY=/leonardo_scratch/large/userinternal/nshukla1/finalhpc/hpc-scaletest/workspace/iPIC3D-CPU-NS/build

date
start_time="$(date -u +%s.%N)"
# Execute command
srun --ntasks-per-node 112 --cpu-bind=cores $BINARY/iPIC3D os-stdin

end_time="$(date -u +%s.%N)"
elapsed="$(bc <<<"$end_time-$start_time")"
echo "Total of $elapsed seconds elapsed for process"
date
