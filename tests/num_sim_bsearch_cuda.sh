#!/bin/bash

source $HOME/src/set_intersection/tests/files.sh

CUDA=$HOME/src/set_intersection/CUDA
#CMD_P=cuda_bsearch
#CMD_P=num_sim_bsearch_gm_cuda
#CMD_P=count_bsearch_gm_cuda
CMD_CS=num_sim_bsearch_sm_cuda
CMD_CG=num_sim_bsearch_gm_cuda
CMD_BG=num_sim_brute_force_gm_cuda

#SEQ=$HOME/src/set_intersection/sequential
#CMD_S=num_sim_scan_seq

#echo "$SEQ/$CMD_S $U $A $B"
#$SEQ/$CMD_S $U $A $B 10000

#echo "$CMD_CGC 1"
#$CUDA/$CMD_CG $U $A $B 10000 1 1024

#echo "$CMD_CGC 1"
#$CUDA/$CMD_CS $U $A $B 10000 1 1024

#U=$HOME/scratch/set_i/hg19.bed

#BED_DIR=$HOME/scratch/set_i/exp

#SEQ=$HOME/src/set_intersection/sequential
#COUNT=num_sim_scan_seq

#echo "$COUNT $i"
#
#for BED_A in $TEST_FILES
#do
#	for BED_B in $TEST_FILES
#	do
#		if [ "$BED_A" != "$BED_B" ]
#		then
#			$CUDA/$CMD_CG $U $BED_A $BED_B 10000 1 1024
#			$CUDA/$CMD_BG $U $BED_A $BED_B 10000 1 1024
#			#echo "$CUDA/$CMD_CG $U $BED_A $BED_B 10000 1 1024"
#		fi
#	done
#done
$CUDA/$CMD_CG $U $A $B 10000 1 1024
$CUDA/$CMD_BG $U $A $B 10000 1 1024
