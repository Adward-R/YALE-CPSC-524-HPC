#!/bin/bash

module load Langs/Intel/15 MPI/OpenMPI/1.8.6-intel15
make clean
make TARGETS=nbody0 TARGETOBJECTS=nbody0.o
make TARGETS=nbody1 TARGETOBJECTS=nbody1.o
make TARGETS=nbody2 TARGETOBJECTS=nbody2.o
make TARGETS=nbody3 TARGETOBJECTS=nbody3.o
make TARGETS=nbody4 TARGETOBJECTS=nbody4.o
