#!/bin/bash

source $HOME/src/set_intersection/tests/files.sh

#CMD_P=cuda_bsearch
#CMD_P=num_sim_bsearch_gm_cuda
#CMD_P=count_bsearch_gm_cuda
CMD_CS=num_sim_bsearch_sm_cuda
CMD_CG=num_sim_bsearch_gm_cuda
CMD_BG=num_sim_brute_force_gm_cuda
CMD_SS=num_sim_scan_seq

for BED_A in $TEST_FILES
do
	for BED_B in $TEST_FILES
	do
		if [ "$BED_A" != "$BED_B" ]
		then
			$CUDA_BASE/$CMD_BG $U $BED_A $BED_B 10000 1 1024
			$SEQ_BASE/$CMD_SS $U $BED_A $BED_B 10000
		fi
	done
done
