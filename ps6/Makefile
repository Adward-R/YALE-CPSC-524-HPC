# This Makefile assumes the following module files are loaded:
#
# Langs/Intel/15
# GPU/Cuda/8.0
#
# This Makefile will only work if executed on a GPU node.
#

CUDAPATH = /home/apps/fas/GPU/Cuda/8.0

NVCC = $(CUDAPATH)/bin/nvcc -Wno-deprecated-gpu-targets 

NVCCFLAGS = -I$(CUDAPATH)/include -O3

LFLAGS = -L$(CUDAPATH)/lib64 -lcuda -lcudart -lm

# Compiler-specific flags (by default, we always use sm_20)
GENCODE_SM20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\"
GENCODE = $(GENCODE_SM20)

.SUFFIXES : .cu .ptx

BINARIES = task1 task2

task1: task1.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<
task2: task2.o
	$(NVCC) $(GENCODE) $(LFLAGS) -o $@ $<


.cu.o:
	$(NVCC) $(GENCODE) $(NVCCFLAGS) -o $@ -c $<

clean:	
	rm -f *.o $(BINARIES)
