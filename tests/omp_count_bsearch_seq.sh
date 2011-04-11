#!/bin/bash

source $HOME/src/set_intersection/tests/files.sh

SEQ=$HOME/src/set_intersection/sequential
COUNT_S=count_bsearch_seq

OMP=$HOME/src/set_intersection/omp
COUNT=count_bsearch_omp


#$OMP/$COUNT $U $A $B


#for i in 3 4 5 6
#for i in 1 2 3 4 5 6
#do
	#echo "$OMP/$COUNT $U $A $B 100 $i"
$OMP/$COUNT $U $A $B 1000 $i 
#done
