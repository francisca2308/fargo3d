#!/bin/bash
#SBATCH --job-name=test_name
#SBATCH --partition=mapu
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --output=/home/fangulo/myforkfargo3d/fargo3d/test.out
#SBATCH --error=/home/fangulo/myforkfargo3d/fargo3d/test.err
module load openmpi/4.1.5
cd /home/fangulo/myforkfargo3d/fargo3d
mpirun -np 4 ./fargo3d setups/dustydisk/dustydisk.par
