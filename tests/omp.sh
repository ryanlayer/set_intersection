#!/bin/bash

#U=/bigtemp/rl6sf/interval/hg18/hg19.bed
#A=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k27me3.broadPeak.bed
#B=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k4me3.broadPeak.bed
#B=$A
#A=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksHsmmH3k4me1.broadPeak.bed
#B=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksH1hescH3k4me1.broadPeak.bed

U=/bigtemp/rl6sf/interval/small/U.bed
A=/bigtemp/rl6sf/interval/small/A.bed
B=/bigtemp/rl6sf/interval/small/B.bed

OMP=$HOME/src/set_intersection/omp
#CMD_P=cuda_bsearch
#CMD_P=num_sim_bsearch_gm_cuda
#CMD_P=count_bsearch_gm_cuda
COUNT=count_bsearch_omp

for i in 6
do
	echo "$COUNT $i"
	$OMP/$COUNT $U $A $B $i 
done
