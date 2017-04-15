#!/bin/bash
#PBS -l procs=1,tpn=1,mem=34gb,walltime=30:00
#PBS -N serial
#PBS -r n
#PBS -q fas_devel
#PBS -j oe

module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15
pwd
cd $PBS_O_WORKDIR
pwd
make clean
make

mpiexec -n 1 ./nbody2

# mpiexec -n 8 ./fserial < ./data/testdata1 > ./testdata1_f.out
# mpiexec -n 8 ./fserial < ./data/testdata2 > ./testdata2_f.out

# mpiexec -n 1 ./serial < ./data/actualdata1 > ./actualdata1_s.out
# mpiexec -n 1 ./serial < ./data/actualdata2 > ./actualdata2_s.out
# mpiexec -n 1 ./serial < ./data/actualdata3 > ./actualdata3_s.out
# mpiexec -n 1 ./serial < ./data/actualdata4 > ./actualdata4_s.out