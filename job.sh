#!/bin/bash
#SBATCH --job-name=sim_2048_2048
#SBATCH --partition=mapu
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --output=/home/fangulo/myforkfargo3d/fargo3d/sim_2048_2048.out
#SBATCH --error=/home/fangulo/myforkfargo3d/fargo3d/sim_2048_2048.err
module load openmpi/4.1.5
cd /home/fangulo/myforkfargo3d/fargo3d
mpirun -np 32 ./fargo3d setups/dustydisk/dustydisk.par
