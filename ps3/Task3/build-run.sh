#!/bin/bash
#PBS -l procs=8,tpn=8,mem=34gb
#PBS -l walltime=30:00
#PBS -N Task3
#PBS -r n
#PBS -j oe
#PBS -q fas_devel

module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15
pwd
cd $PBS_O_WORKDIR
pwd
cat $PBS_NODEFILE
make clean
make task3
time mpiexec -n 4 ./task3 4000
time mpiexec -n 8 ./task3 8000
