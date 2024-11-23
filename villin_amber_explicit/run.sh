#!/bin/bash
#SBATCH -p all
#SBATCH -w n5
#SBATCH -J Villin_explicit
#SBATCH -n 1  #num of total mpi processes
#SBATCH -c 2  #num of threads per mpi processes
#SBATCH --mail-type=ALL
#SBATCH -o run.log

# set GPU ID
export CUDA_VISIBLE_DEVICES=3

# execution
#python3 1_build.py
python3 2_equilibration.py
python3 3_production.py

