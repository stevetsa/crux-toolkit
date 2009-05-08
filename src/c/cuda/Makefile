CUDA_CC = nvcc
CC=gcc
CUDA_COPTS = --host-compilation 'c' -DCUDA_NVCC --device-emulation
CUDA_DIR=/usr/local/cuda
#CUDA_DIR=/net/gs/vol3/software/modules-sw/cuda/2.1/Linux/RHEL5/i686/


all: crux_cuda_test


crux_cuda_test: crux_cuda_test.o crux_cuda.cu.o
	$(CC) crux_cuda_test.o crux_cuda.cu.o -L$(CUDA_DIR)/lib -lcuda -lcudart -o crux_cuda_test


crux_cuda_test.o:crux_cuda_test.c
	$(CC) -c crux_cuda_test.c -o crux_cuda_test.o

crux_cuda.o:crux_cuda.c crux_cuda.h
	$(CC) -c crux_cuda.c -o crux_cuda.o

crux_cuda.cu.o:crux_cuda.cu
	$(CUDA_CC) $(CUDA_COPTS) -c crux_cuda.cu -o crux_cuda.cu.o


clean:
	rm *.o crux_cuda_test *~ crux_cuda.linkinfo
