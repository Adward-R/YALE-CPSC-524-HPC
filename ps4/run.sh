#!/bin/bash
#PBS -l procs=8,tpn=2,mem=34gb,walltime=30:00
#PBS -N fserial
#PBS -r n
#PBS -q fas_devel
#PBS -j oe

module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15
pwd
cd $PBS_O_WORKDIR
pwd
make clean
make

# mpiexec -n 8 ./fserial < ./data/testdata1 > ./testdata1_f.out
# mpiexec -n 8 ./fserial < ./data/testdata2 > ./testdata2_f.out

mpiexec -n 8 ./fserial < ./data/actualdata1 > ./actualdata1_f.out
mpiexec -n 8 ./fserial < ./data/actualdata2 > ./actualdata2_f.out
mpiexec -n 8 ./fserial < ./data/actualdata3 > ./actualdata3_f.out
mpiexec -n 8 ./fserial < ./data/actualdata4 > ./actualdata4_f.out