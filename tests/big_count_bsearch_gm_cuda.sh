#!/bin/bash

source $HOME/src/set_intersection/tests/files.sh

BASE=/localtmp/rl6sf

EXON=$BASE/exon_hg19.bed
K=$BASE/NA12878.chr1.1K.bed
TEN_K=$BASE/NA12878.chr1.10K.bed
HUN_K=$BASE/NA12878.chr1.100K.bed
M=$BASE/NA12878.chr1.1M.bed
TEN_M=$BASE/NA12878.chr1.10M.bed
HUN_M=$BASE/NA12878.chr1.100M.bed

#$COUNT_BSEARCH_GM_CUDA $U $TEN_KIL_A $TEN_KIL_B 1 1024 0
#$COUNT_BSEARCH_GM_CUDA $U $HUND_KIL_A $HUND_KIL_B 1 1024 0
#$COUNT_BSEARCH_GM_CUDA $U $MIL_A $MIL_B 1 1024 1
#$COUNT_BSEARCH_GM_CUDA $U $TEN_MIL_A $TEN_MIL_B 1 1024 0

#$BIG_COUNT_BSEARCH_GM_CUDA $U $MIL_A $MIL_B 1 1024 1000
#$BIG_COUNT_BSEARCH_GM_CUDA $U $MIL_A $MIL_B 1 1024 10000
#$BIG_COUNT_BSEARCH_GM_CUDA $U $EXON $NA 1 1024 1000000 1

DEV=1

$BIG_COUNT_BSEARCH_GM_CUDA $U $EXON $K 1 1024 1000 $DEV

$BIG_COUNT_BSEARCH_GM_CUDA $U $EXON $HUN_K 1 1024 10000 $DEV

$BIG_COUNT_BSEARCH_GM_CUDA $U $EXON $M 1 1024 100000 $DEV

$BIG_COUNT_BSEARCH_GM_CUDA $U $EXON $TEN_M 1 1024 1000000 $DEV

$BIG_COUNT_BSEARCH_GM_CUDA $U $EXON $HUN_M 1 1024 10000000 $DEV