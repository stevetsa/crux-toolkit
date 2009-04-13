#ifndef _CRUX_CUDA_CU_
#define _CRUX_CUDA_CU_

#include "crux_cuda.cu.h"

#include <stdio.h>

#define MAX_XCORR_OFFSET 50
#define NUM_THREADS_PER_BLOCK 256
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>


/*
__global__ static void reduction(float*a, int stride) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int n = stride * 2;
  if (i * n < N) {
    int index = i * n;
    a[index] = max(a[index], a[index] + stride);
  }

}
*/

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



void cuda_normalize_each_region(float* h_values, float* max_per_region, int n, int num_regions, int region_selector) {
  float *d_values;
  float *d_max_per_region;

  cudaError error = cudaMalloc((void**)&d_values, sizeof(float)*n);
  if (error != cudaSuccess)
    {
      printf("Error allocating d_values:%i",error);
    }
  
  error = cudaMalloc((void**)&d_max_per_region, sizeof(float)* num_regions);
  
  error = cudaMemcpy(d_values, h_values, sizeof(float)*n, cudaMemcpyHostToDevice);
  error = cudaMemcpy(d_max_per_region, max_per_region, sizeof(float)*n, cudaMemcpyHostToDevice);

  int num_blocks = n / NUM_THREADS_PER_BLOCK;

  d_normalize_each_region<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_values, d_max_per_region, n, num_regions, region_selector);
  
  error = cudaGetLastError();

  cudaMemcpy(h_values, d_values, sizeof(float)*n, cudaMemcpyDeviceToHost);

}


/***********************************************
 *cross_correlation_obs(): Intermediate function
 *to learn about kernel functions, will do the
 *cross_correlation part of the observed spectrum
 *on the gpu.  Much can be done to optimize this, i.e
 *the next step would be to do all of the observed data
 *processing on the card via a specialized kernel
 *function.
 *
 *The steps to take would be
 *1) Load raw data onto card.
 *2) Sqrt data
 *3) find max per region
 *4) normalize to 50 for each region.
 *5) do cross correlation algorithm.

 ***********************************************/
void cross_correlation_obs(float* h_values, float* h_ans, int n, int max_offset) {
  float *d_values;
  float *d_ans;

  cudaError error = cudaMalloc((void**)&d_values, sizeof(float)*n);
  if (error != cudaSuccess) 
    printf("Error allocating d_values:%i",error);

  error = cudaMalloc((void**)&d_ans, sizeof(float)*n);
  if (error != cudaSuccess) 
    printf("Error allocating d_ans:%i",error);

  
  error = cudaMemcpy(d_values, h_values, sizeof(float)*n, cudaMemcpyHostToDevice);
  if (error != cudaSuccess)
    printf("Error copying values:%i",error);

  int num_blocks = n / NUM_THREADS_PER_BLOCK;

  d_cross_correlation_obs<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_values, d_ans, n, max_offset);

  //check for errors.
  error = cudaGetLastError();

  if (error != cudaSuccess)
    printf ("Kernel error:%i\n",error);
  
  //copy answer back.
  cudaMemcpy(h_ans, d_ans, sizeof(float)*n, cudaMemcpyDeviceToHost);

}
 






#endif
