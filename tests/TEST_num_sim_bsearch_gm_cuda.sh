#!/bin/bash

source $HOME/src/set_intersection/tests/files.sh


I=10

for i in {1..10}
do
$SIM_BSEARCH_GM_CUDA $U $KIL_A $KIL_B 1000 1 1024 0
$SIM_BSEARCH_GM_CUDA $U $TEN_KIL_A $TEN_KIL_B 1000 1 1024 0
$SIM_BSEARCH_GM_CUDA $U $HUND_KIL_A $HUND_KIL_B 1000 1 1024 0
$SIM_BSEARCH_GM_CUDA $U $MIL_A $MIL_B 1000 1 1024 1
$SIM_BSEARCH_GM_CUDA $U $TEN_MIL_A $TEN_MIL_B 1000 1 1024 0
done