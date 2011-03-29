#!/bin/bash

U=/bigtemp/rl6sf/interval/hg18/hg19.bed
A=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k27me3.broadPeak.bed
B=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k4me3.broadPeak.bed
#B=$A
#A=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksHsmmH3k4me1.broadPeak.bed
#B=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksH1hescH3k4me1.broadPeak.bed

#U=/bigtemp/rl6sf/interval/small/U.bed
#A=/bigtemp/rl6sf/interval/small/A.bed
#B=/bigtemp/rl6sf/interval/small/B.bed

CUDA=$HOME/src/set_intersection/CUDA
#CMD_P=cuda_bsearch
#CMD_P=num_sim_bsearch_gm_cuda
#CMD_P=count_bsearch_gm_cuda
CMD_CGE=enumerate_bsearch_gm_cuda
CMD_CGC=count_bsearch_gm_cuda
CMD_CSC=count_bsearch_sm_cuda

SEQ=$HOME/src/set_intersection/sequential
CMD_S=$SEQ/enumerate_bsearch_seq

echo "$CMD_CGC 1"
$CUDA/$CMD_CGC $U $A $B 1 1 1024 1

echo "$CMD_CGC 0"
$CUDA/$CMD_CGC $U $A $B 1 1 1024 0

echo "$CMD_CSC 1"
$CUDA/$CMD_CGC $U $A $B 1 1 1024 1

echo "$CMD_CSC 0"
$CUDA/$CMD_CGC $U $A $B 1 1 1024 0
