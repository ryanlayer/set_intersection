#!/bin/bash

source $HOME/src/set_intersection/tests/files.sh

CUDA=$HOME/src/set_intersection/CUDA
#CMD_P=cuda_bsearch
#CMD_P=num_sim_bsearch_gm_cuda
#CMD_P=count_bsearch_gm_cuda
CMD_CGE=enumerate_bsearch_gm_cuda
CMD_CGC=count_bsearch_gm_cuda
CMD_CSC=count_bsearch_sm_cuda

SEQ=$HOME/src/set_intersection/sequential
CMD_S=count_scan_seq

#echo "$SEQ/$CMD_S $U $A $B"
$SEQ/$CMD_S $U $A $B

#echo "$CMD_CGC 1"
$CUDA/$CMD_CGC $U $A $B 1 1 1024 1

#echo "$CMD_CGC 0"
$CUDA/$CMD_CGC $U $A $B 1 1 1024 0

#echo "$CMD_CSC 1"
$CUDA/$CMD_CGC $U $A $B 1 1 1024 1

#echo "$CMD_CSC 0"
$CUDA/$CMD_CGC $U $A $B 1 1 1024 0
