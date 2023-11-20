#!/bin/bash
#SBATCH -p all
#SBATCH -J alanine_dipeptide # job name
#SBATCH -n 1  # num of total mpi processes
#SBATCH -c 1  # num of threads per mpi processes
#SBATCH -o run.log

# set GPU ID
#export CUDA_VISIBLE_DEVICES=0

python sim.py alanine-dipeptide-nowater.pdb GB99dms.xml traj
#周期境界条件を用いる場合
#python sim.py alanine-dipeptide-nowater.pdb GB99dms.xml traj CutoffPeriodic