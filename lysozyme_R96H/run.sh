#!/bin/bash
#SBATCH -p all
#SBATCH -J R96H
#SBATCH -n 1  #num of total mpi processes
#SBATCH -c 4  #num of threads per mpi processes
#SBATCH -o run.log

# set GPU ID
#export CUDA_VISIBLE_DEVICES=0

# execution
#python3 1_build.py
python3 2_equilibration.py
python3 3_production.py

