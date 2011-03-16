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

#CMD_P=../../bin/linux/release/cuda_bsearch
CMD_P=cuda_bsearch
#CMD_P=../../bin/linux/release/cuda_num_sim
CMD_S1=$HOME/src/set_intersection/sequential/bsearch_seq
#CMD_S2=$HOME/src/tools/bed/cluster.pl
#CMD_S3=$HOME/src/set_intersection/sequential/scan_seq
#CMD_S4=$HOME/src/set_intersection/sequential/brute_force

echo $CMD_P
$CMD_P $U $A $B $1
#$CMD_S1 $U $A $B $1
#$CMD_P $U $A $B 10| awk '{sum+=$1} END {print sum}'

#echo $CMD_S1
#$CMD_S1 $U $A $B 10

#echo $CMD_S3
#$CMD_S3 $U $A $B

#echo $CMD_S4
#$CMD_S4 $U $A $B

#A=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksHsmmH3k4me1.broadPeak.bed
#B=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksH1hescH3k4me1.broadPeak.bed
#$CMD_P $U $A $B $1
#$CMD_S1 $U $A $B $1
