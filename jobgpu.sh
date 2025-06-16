#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=test.out
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --nodelist=g001
#SBATCH -p kurruf_gpu
#SBATCH --gres=gpu:A100:1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load openmpi

mpirun -np 1 ./fargo3d setups/dustydisk/dustydisk.par