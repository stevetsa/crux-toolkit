
#ifndef _CRUX_CUDA_CU_
#define _CRUX_CUDA_CU_

#include "crux_cuda.cu.h"

#include <stdio.h>

#define NUM_THREADS_PER_BLOCK 256

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
//#include "carp.h"

void calcMax(int n, float*d_ans);
void calcMaxP2(float* d_in, int n, float* d_ans);
void calcMaxNP2(float* d_in, int n, float* d_ans);


/***************************************************************
 * KERNELS
 ***************************************************************/
__global__ static void reduction(float*a, int stride, int N) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int n = stride * 2;
  if (i * n < N) {
    int index = i * n;
    a[index] = max(a[index], a[index + stride]);
  }

}


__global__ static void d_sqrt(float* invalues) {
  const int idx = threadIdx.x + blockIdx.x*NUM_THREADS_PER_BLOCK;
  float ans = sqrt(invalues[idx]);
  invalues[idx] = ans;
}

__global__ static void d_normalize_each_region(float* invalues, 
					       float* max_per_region, 
					       int max_N, 
					       int num_regions,
					       int region_selector) {

  const int idx = threadIdx.x + blockIdx.x*NUM_THREADS_PER_BLOCK;
  int region = idx / region_selector;
  if (max_per_region[region] == 0) return;
  invalues[idx] = invalues[idx] / max_per_region[region] * 50;
}


__global__ static void d_cross_correlation_obs(float* invalues, float* ans, int max_N, int max_offset) {
  const int idx = threadIdx.x + blockIdx.x*NUM_THREADS_PER_BLOCK;
  
  int cur_idx;
  int min_idx = max(0, idx - max_offset);
  int max_idx = min(max_N-1, idx + max_offset);

  ans[idx] = invalues[idx];
  
  for (cur_idx = min_idx ; cur_idx <= max_idx ; cur_idx++) {
    ans[idx] -= (invalues[cur_idx] / (max_offset * 2.0));
  }
}

/***************************************************************
 * Device static variables
 ***************************************************************/
float* d_temp = NULL;
int d_temp_size = -1;
float* d_values = NULL;
int d_values_size = -1;
float* d_max_per_region = NULL;
int d_max_per_region_size = -1;
float* d_corr = NULL;
int d_corr_size = -1;

float* d_ans = NULL;
int d_ans_size = -1;

/***************************************************************
 * HOST FUNCTIONS AND PROCEDURES
 ***************************************************************/
void cuda_float_malloc(float** a, int* current_size, int size) {
  if (*current_size < size) {
    if (*a == NULL)
      CUDAEXEC(cudaFree(*a),"free(*a)");
    CUDAEXEC(cudaMalloc((void**)a, sizeof(float) * size),"alloc(*a)");
    *current_size = size;
  }
}

void crux_cuda_initialize() {
  printf("crux_cuda_initialize(): start\n");
  cuda_float_malloc(&d_values, &d_values_size, 1024);
  cuda_float_malloc(&d_temp, &d_temp_size, 1024);
  cuda_float_malloc(&d_ans, &d_ans_size, 1);
  cuda_float_malloc(&d_max_per_region, &d_max_per_region_size, 10);
  cuda_float_malloc(&d_corr, &d_corr_size, 1024);
}

void crux_cuda_shutdown() {
  CUDAEXEC(cudaFree(d_values),"free(d_values)");
  CUDAEXEC(cudaFree(d_temp),"free(d_temp)");
  CUDAEXEC(cudaFree(d_ans),"free(d_ans)");
}



int isPowerOfTwo(int n)
{
    return (n) && !(n & (n - 1)); //this checks if the integer n is a power of two or not
}

/*
 *calcMax: so the reduction code works for arrays with
 *elements in powers of 2.
 *So my current solution is to copy the array into
 *an array with a power of two size and run the algorithm
 *I am assuming zero is the smallest value that the intensity
 *array can take (which should be true since we are sqrt
 *the intensties).  Later I should take a look at 
 *improving the efficiency, but there is an example of
 *a sum reduction on the web that is much more efficient
 *than what I current have.  I will have to study that
 *and incorporate it into the code, since this reduction
 *is the rate limiting step.  I wish I knew more about 
 *parallelism and calculating cost and such.
 *****************************************************/
void calcMax(float* d_in, int n, float*d_ans) {

  if (isPowerOfTwo(n)) {
    calcMaxP2(d_in, n, d_ans);
  }
  else {
    calcMaxNP2(d_in, n, d_ans);
  }
}

void doReductionMax(float* d_in, int n, float* d_ans) {
  int i;
  int stride;
  int num_threads;
  int num_blocks;
  cudaError error;

  if (n < NUM_THREADS_PER_BLOCK) {
    num_threads = n;
    num_blocks = 1;
  }
  else {
    num_threads = NUM_THREADS_PER_BLOCK;
    num_blocks = n / num_threads;
    if (n % num_threads !=0) num_blocks++;
  }

  for(stride=1;stride<=n/2;stride*=2) {
    reduction<<<num_blocks, num_threads>>>(d_in, stride, n);
    error = cudaGetLastError();
    if (error != cudaSuccess)
      printf("reduction error %s stride:%i n:%i\n",cudaGetErrorString(error),stride, n);
    cudaThreadSynchronize();
  }

  CUDAEXEC(cudaMemcpy(d_ans, 
		      d_in, 
		      sizeof(float), 
		      cudaMemcpyDeviceToDevice),
	   "d_values -> d_ans");  
}

/******************************************************
 *calcMaxP2: calculates the max in an array that has
 *a power of 2 elements
 ******************************************************/
void calcMaxP2(float*d_in, int n, float* d_ans) {
  cuda_float_malloc(&d_temp, &d_temp_size, n);
  //printf("n:%d num_threads:%d num_blocks:%d\n",n,num_threads, num_blocks);

  CUDAEXEC(cudaMemcpy(d_temp, d_in, n*sizeof(float), cudaMemcpyDeviceToDevice),
	   "d_temp <- d_value");
  
  doReductionMax(d_temp, n, d_ans);
}


/*********************************************************
 *calcMaxNP2: calculates the max in an array that does not
 *have power of 2 elements. It achieves this by copying to
 *a temporary array that is a power of 2.
 *********************************************************/
void calcMaxNP2(float* d_in, int n, float* d_ans){
  int n2 = (int)pow(2,ceil(log2((float)n)));

  cuda_float_malloc(&d_temp, &d_temp_size, n2);
  CUDAEXEC(cudaMemset(d_temp,0,n2*sizeof(float)),"memset(d_temp,0)");
  CUDAEXEC(cudaMemcpy(d_temp,
		      d_in,
		      n*sizeof(float),
		      cudaMemcpyDeviceToDevice),
	   "calcMaxNP2: d_temp <- d_values");
  doReductionMax(d_temp, n2, d_ans);
}

void calcMax2(float* d_in, int start, int end, float* d_ans) {
  if (end-start < 1) {
    cudaMemset(d_ans, 0, sizeof(float));
  }
  else {
    calcMax(d_in+start,end-start,d_ans);
  }
}


void cuda_sqrt_max_normalize_and_cc(float* h_values, int n, int num_regions, 
				    int region_selector, int max_offset) {
  int size_n = n * sizeof(float);

  //printf("allocating memory %i %i\n",size_n, size_reg);

  cuda_float_malloc(&d_values, &d_values_size, n);

  CUDAEXEC(cudaMemcpy(d_values, h_values, 
		      size_n, cudaMemcpyHostToDevice),
	   "h_values -> d_values");

  d_cuda_sqrt_max_normalize_and_cc(d_values, n, num_regions, region_selector, max_offset);

  CUDAEXEC(cudaMemcpy(h_values, d_corr, size_n, cudaMemcpyDeviceToHost),"d_corr -> h_values");

}


void d_cuda_sqrt_max_normalize_and_cc2(float* d_in, float* d_out, int n, int num_regions,
				       int region_selector, int max_offset) {
  cudaError error;
  int region;
  cuda_float_malloc(&d_max_per_region, &d_max_per_region_size, num_regions);
 
  int num_blocks = n / NUM_THREADS_PER_BLOCK; 

  //printf("Executing d_sqrt_and_max_region\n");
  d_sqrt<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_in); 
  error = cudaGetLastError(); 
  if (error != cudaSuccess) 
    printf("d_sqrt: error %s\n",cudaGetErrorString(error)); 

  for(region=0;region<num_regions;region++) {
    int start = region_selector * region;
    int end = min(n, start + region_selector);
    calcMax2(d_in, start, end, d_max_per_region+region);
  }

  d_normalize_each_region<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_in, d_max_per_region,  
                                                                 n, num_regions, region_selector); 
  error = cudaGetLastError(); 
  
  if (error != cudaSuccess) {
    printf("d_normalize_each_region: error %s\n",cudaGetErrorString(error));
    printf("n: %i nr: %i s: %i nb: %i ntb: %i\n",n,num_regions,region_selector, 
	   num_blocks, NUM_THREADS_PER_BLOCK);
  }
  d_cross_correlation_obs<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_in, d_out, n, max_offset); 
  error = cudaGetLastError(); 

  if (error != cudaSuccess)
    printf("d_cross_correlation_obs: error %s\n", cudaGetErrorString(error));

}

void d_cuda_sqrt_max_normalize_and_cc(float* d_in, int n, int num_regions, 
				      int region_selector, int max_offset) {
  cuda_float_malloc(&d_corr, &d_corr_size, n);
  d_cuda_sqrt_max_normalize_and_cc2(d_in, d_corr, n, num_regions, region_selector, max_offset);
}




#endif
