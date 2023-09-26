#!/bin/bash
#SBATCH -p all
#SBATCH -J prod
#SBATCH -n 1  #num of total mpi processes
#SBATCH -c 2  #num of threads per mpi processes
#SBATCH -o run.log

# set GPU ID
#export CUDA_VISIBLE_DEVICES=0

# execution
/opt/namd/namd3 +p 2 run.inp

