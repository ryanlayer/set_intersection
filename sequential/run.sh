#!/bin/bash

#U=/bigtemp/rl6sf/interval/hg18/hg18.bed
#A=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k27me3.broadPeak.bed
#B=/bigtemp/rl6sf/interval/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k4me3.broadPeak.bed
#B=$A

U=/bigtemp/rl6sf/interval/small/U.bed
A=/bigtemp/rl6sf/interval/small/A.bed
B=/bigtemp/rl6sf/interval/small/B.bed

$1 $U $A $B 10
