#!/bin/bash
#SBATCH -p all
#SBATCH -w n5
#SBATCH -J villin_implicit # job name
#SBATCH -n 1  # num of total mpi processes
#SBATCH -c 1  # num of threads per mpi processes
#SBATCH --mail-type=all
#SBATCH -o run.log

# set GPU ID
export CUDA_VISIBLE_DEVICES=3

#python sim.py alanine-dipeptide-nowater.pdb amber traj
#周期境界条件を用いる場合
python sim.py system.pdb amber traj CutoffPeriodic
