#!/bin/bash

source $HOME/src/set_intersection/tests/files.sh

REPS=1000

function run {

for i in 1 2 3 4 5 6 7 8
do
	echo $i
	$SIM_BSEARCH_OMP $1 $2 $3 $4 $i
done
}


run $U $KIL_A $KIL_B $REPS 
run $U $TEN_KIL_A $TEN_KIL_B $REPS 
run $U $HUND_KIL_A $HUND_KIL_B $REPS 
run $U $MIL_A $MIL_B $REPS 
run $U $TEN_MIL_A $TEN_MIL_B $REPS 
