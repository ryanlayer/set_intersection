#!/bin/bash

source $HOME/src/set_intersection/tests/files.sh

SEQ=$HOME/src/set_intersection/sequential
COUNT_S=count_bsearch_seq

MPI=$HOME/src/set_intersection/mpi
COUNT=num_sim_bsearch_mpi


#$OMP/$COUNT $U $A $B


#for i in 3 4 5 6
for i in 1 2 3 4 5 6
do
	#echo "$OMP/$COUNT $U $A $B 100 $i"
	time mpiexec -n $i $MPI/$COUNT $U $A $B 1000
done
