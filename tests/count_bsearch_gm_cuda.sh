#!/bin/bash

source $HOME/src/set_intersection/tests/files.sh

$COUNT_BSEARCH_GM_CUDA $U $KIL_A $KIL_B 1 1024 0
$COUNT_BSEARCH_GM_CUDA $U $TEN_KIL_A $TEN_KIL_B 1 1024 0
$COUNT_BSEARCH_GM_CUDA $U $HUND_KIL_A $HUND_KIL_B 1 1024 0
$COUNT_BSEARCH_GM_CUDA $U $MIL_A $MIL_B 1 1024 1
$COUNT_BSEARCH_GM_CUDA $U $TEN_MIL_A $TEN_MIL_B 1 1024 0
