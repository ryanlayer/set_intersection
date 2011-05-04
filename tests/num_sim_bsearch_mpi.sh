#!/bin/bash
#PBS -l select=30:mem=2gb
#PBS -l walltime=10:00:00
# PBS -j oe
#PBS -m ae

cd $PBS_O_WORKDIR
source $HOME/src/set_intersection/tests/files.sh

#time mpiexec -n 5 $SIM_BSEARCH_MPI $U $HUND_KIL_A $HUND_KIL_B 1000
mpiexec -n 30 $SIM_BSEARCH_MPI $U $MIL_A $MIL_B 1000
