#TESTS=1 5 10 15 20 25
TESTS="5 20"
for i in $TESTS
do
echo "#!/bin/bash" > run.$i.pbs
echo "#PBS -l select=$i:mem=2gb" >> run.$i.pbs
echo "#PBS -l walltime=10:00:00" >> run.$i.pbs
echo "#PBS -m ae" >> run.$i.pbs
echo "#PBS -o num_sim_bsearch_mpi.run.n-$i.out" >> run.$i.pbs
echo "cd \$PBS_O_WORKDIR" >> run.$i.pbs
echo "source \$HOME/src/set_intersection/tests/files.sh" >> run.$i.pbs
echo "mpiexec -n $i \$SIM_BSEARCH_MPI \$U \$HUND_KIL_A \$HUND_KIL_B 1000" >> run.$i.pbs
#qsub  run.$i.pbs
#rm  run.$i.pbs
done
