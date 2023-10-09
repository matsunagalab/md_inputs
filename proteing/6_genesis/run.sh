#!/bin/bash

export PATH=/opt/genesis-2.1.1/bin:$PATH
export OMP_NUM_THREADS=1

mpirun -np 8 spdyn INP

