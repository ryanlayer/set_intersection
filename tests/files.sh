#!/bin/bash

CUDA_BASE=$HOME/src/set_intersection/CUDA
SEQ_BASE=$HOME/src/set_intersection/sequential
OMP_BASE=$HOME/src/set_intersection/omp
MPI_BASE=$HOME/src/set_intersection/mpi

#BASE=/home/rl6sf/data/intr
BASE=/bigtemp/rl6sf/interval/

U=$BASE/hg18/hg19.bed
A=$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k27me3.broadPeak.bed
B=$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k4me3.broadPeak.bed
BIGA=$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep1Gm12865Ctcf.broadPeak.bed
BIGB=$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep1Gm12865Ctcf.broadPeak.bed
SIMPLE_REPEATS=$BASE/hg18/simple_repeats.bed
#B=$A
#A=$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksHsmmH3k4me1.broadPeak.bed
#B=$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksH1hescH3k4me1.broadPeak.bed

#U=$BASE/small/U.bed
#A=$BASE/small/A.bed
#B=$BASE/small/B.bed

TEST_FILES="$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878Ctcf.broadPeak.bed
$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k27ac.broadPeak.bed
$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k27me3.broadPeak.bed
$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k36me3.broadPeak.bed
$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k4me1.broadPeak.bed
$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k4me2.broadPeak.bed
$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k4me3.broadPeak.bed
$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k9ac.broadPeak.bed
$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H4k20me1.broadPeak.bed
$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep1Gm12878Ctcf.broadPeak.bed
$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep1Gm12878H3k27me3.broadPeak.bed
$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep1Gm12878H3k36me3.broadPeak.bed
$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep1Gm12878H3k4me3.broadPeak.bed
$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep2Gm12878Ctcf.broadPeak.bed
$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep2Gm12878H3k27me3.broadPeak.bed
$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep2Gm12878H3k36me3.broadPeak.bed
$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep2Gm12878H3k4me3.broadPeak.bed"
