#!/bin/bash

export PATH=/opt/genesis-2.1.0/bin:$PATH
export OMP_NUM_THREADS=1

mpirun -np 8 spdyn INP1
mpirun -np 8 spdyn INP2
mpirun -np 8 spdyn INP3

