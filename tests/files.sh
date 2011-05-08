#!/bin/bash

CUDA_BASE=$HOME/src/set_intersection/CUDA
SEQ_BASE=$HOME/src/set_intersection/sequential
OMP_BASE=$HOME/src/set_intersection/omp
MPI_BASE=$HOME/src/set_intersection/mpi

COUNT_SWEEP_SEQ=$SEQ_BASE/count_sweep_seq
COUNT_BRUTE_SEQ=$SEQ_BASE/count_brute_force_seq
COUNT_BSEARCH_SEQ=$SEQ_BASE/count_bsearch_seq
BIG_COUNT_BSEARCH_SEQ=$SEQ_BASE/big_count_bsearch_seq
SIM_BSEARCH_SEQ=$SEQ_BASE/num_sim_bsearch_seq

COUNT_BSEARCH_GM_CUDA=$CUDA_BASE/count_bsearch_gm_cuda
SIM_BSEARCH_GM_CUDA=$CUDA_BASE/num_sim_bsearch_gm_cuda

COUNT_BSEARCH_OMP=$OMP_BASE/count_bsearch_omp
SIM_BSEARCH_OMP=$OMP_BASE/num_sim_bsearch_omp

SIM_BSEARCH_MPI=$MPI_BASE/num_sim_bsearch_mpi

#BASE=/home/rl6sf/data/intr
#BASE=$HOME/data/intervals/
#BASE=/bigtemp/rl6sf/interval/
BASE=/home/rl6sf/data/interval

U=$BASE/hg18/hg19.bed
A=$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k27me3.broadPeak.bed
B=$BASE/hg18/broad/wgEncodeBroadChipSeqPeaksGm12878H3k4me3.broadPeak.bed
A_a=$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep1Gm12878H3k27me3.broadPeak.bed
A_b=$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep1Gm12878H3k27me3.broadPeak.bed
BIGA=$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep1Gm12865Ctcf.broadPeak.bed
BIGB=$BASE/hg18/UW/wgEncodeUwChIPSeqHotspotsRep1Gm12865Ctcf.broadPeak.bed
SIMPLE_REPEATS=$BASE/hg18/simple_repeats.bed
TEN_MIL_A=$BASE/hg18/c.aln.10M.fix.bed
TEN_MIL_B=$BASE/hg18/c.aln.10M.shuff.fix.bed
MIL_A=$BASE/hg18/c.aln.1M.fix.bed
MIL_B=$BASE/hg18/c.aln.1M.shuff.fix.bed
HUND_KIL_A=$BASE/hg18/c.aln.100K.fix.bed
HUND_KIL_B=$BASE/hg18/c.aln.100K.shuff.fix.bed
TEN_KIL_A=$BASE/hg18/c.aln.10K.fix.bed
TEN_KIL_B=$BASE/hg18/c.aln.10K.shuff.fix.bed
KIL_A=$BASE/hg18/c.aln.1K.fix.bed
KIL_B=$BASE/hg18/c.aln.1K.shuff.fix.bed
HUND_A=$BASE/hg18/c.aln.100.fix.bed
HUND_B=$BASE/hg18/c.aln.100.shuff.fix.bed

TINY_U=$BASE/small/oU.bed
TINY_A=$BASE/small/oA.bed
TINY_B=$BASE/small/oB.bed

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
