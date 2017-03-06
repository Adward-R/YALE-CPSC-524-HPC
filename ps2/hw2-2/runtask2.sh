#!/bin/bash
#PBS -l procs=4,tpn=2,mem=68gb
#PBS -l walltime=2:00
#PBS -N Task1
#PBS -r n
#PBS -j oe
#PBS -q cpsc424

# Load necessary module files
module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15

# Print initial working directory
pwd

# Change to submission directory
cd $PBS_O_WORKDIR

# Now print the current working directory
pwd

# Print the node list
cat $PBS_NODEFILE

# Run the program 3 times
time mpiexec -n 4 task2
time mpiexec -n 4 task2
time mpiexec -n 4 task2
