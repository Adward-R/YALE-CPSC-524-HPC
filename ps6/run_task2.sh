#!/bin/bash
#PBS -l procs=5,tpn=5,mem=40gb,walltime=30:00
#PBS -N massive
#PBS -r n
#PBS -j oe
#PBS -q cpsc424gpu

module load Langs/Intel/15 GPU/Cuda/8.0
pwd
cd $PBS_O_WORKDIR
pwd
device=2
./task2 1024 1024 1024 $device
./task2 8192 8192 8192 $device
./task2 1024 8192 1024 $device
./task2 8192 1024 8192 $device
./task2 8192 8192 1024 $device
# ./task1 8192 8192 8192 2
# ./task1 1024 1024 1024 2
