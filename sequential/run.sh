#!/bin/bash

#U=/bigtemp/rl6sf/interval/hg18/hg18.bed
#A=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k27me3.broadPeak.bed
#B=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k4me3.broadPeak.bed
#B=$A

U=/bigtemp/rl6sf/interval/hg18/hg19.bed
#A=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k27me3.broadPeak.bed
#B=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k4me3.broadPeak.bed
#B=$A
A=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksHsmmH3k4me1.broadPeak.bed
B=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksH1hescH3k4me1.broadPeak.bed


#U=/bigtemp/rl6sf/interval/small/U.bed
#A=/bigtemp/rl6sf/interval/small/A.bed
#B=/bigtemp/rl6sf/interval/small/B.bed


echo $1 $U $A $B 10
$1 $U $A $B 10

A=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k27me3.broadPeak.bed
B=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k4me3.broadPeak.bed

echo $1 $U $A $B 10
$1 $U $A $B 10
