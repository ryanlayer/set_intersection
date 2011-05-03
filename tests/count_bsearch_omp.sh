#!/bin/bash

source $HOME/src/set_intersection/tests/files.sh


#for i in 1 2 3 4 5 6 7 8
#do
#$COUNT_BSEARCH_OMP $U $KIL_A $KIL_B $i
#done
#echo
#
#for i in 1 2 3 4 5 6 7 8
#do
#$COUNT_BSEARCH_OMP $U $TEN_KIL_A $TEN_KIL_B 1
#done
#
#for i in 1 2 3 4 5 6 7 8
#do
#$COUNT_BSEARCH_OMP $U $HUND_KIL_A $HUND_KIL_B 1
#done
#$COUNT_BSEARCH_OMP $U $HUND_KIL_A $HUND_KIL_B 
#
#for i in 1 2 3 4 5 6 7 8
#do
#$COUNT_BSEARCH_OMP $U $MIL_A $MIL_B 5
#done
#
for i in 3 4 5 6 7 8
do
echo $i
$COUNT_BSEARCH_OMP $U $TEN_MIL_A $TEN_MIL_B $i
done
