#!/bin/bash
# qsub -I -l procs=5,tpn=5,mem=40gb,walltime=1:00:00 -q cpsc424gpu
module load Langs/Intel/15 GPU/Cuda/8.0
cd ~/YALE-CPSC-524-HPC/ps6/
~ahs3/bin/gpuget
~ahs3/bin/gpulist
