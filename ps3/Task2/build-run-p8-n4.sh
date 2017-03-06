#!/bin/bash
#PBS -l procs=8,tpn=2,mem=34gb
#PBS -l walltime=30:00
#PBS -N Task2
#PBS -r n
#PBS -j oe
#PBS -q cpsc424

module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15
pwd
cd $PBS_O_WORKDIR
pwd
cat $PBS_NODEFILE
make clean
make task2
time mpiexec -n 8 ./task2 8000
