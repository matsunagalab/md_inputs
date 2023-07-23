#!/bin/bash
#PJM -L "node=32"
#PJM -L "rscgrp=small"
#PJM -L "elapse=10:00:00"
#PJM --mpi "max-proc-per-node=4"
#PJM –x PJM_LLIO_GFSCACHE=/vol0004
#PJM -g hp230209
#PJM -j
#PJM -S

export OMP_NUM_THREADS=12
export PLE_MPI_STD_EMPTYFILE=off

bindir=/vol0403/data/hp230209/matsunaga/genesis-2.1.0/bin
mpiexec --of-proc run.out ${bindir}/spdyn INP

