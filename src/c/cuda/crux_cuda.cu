#ifndef _CRUX_CUDA_CU_
#define _CRUX_CUDA_CU_

#include "crux_cuda.cu.h"

#include <stdio.h>

#define MAX_XCORR_OFFSET 50
#define NUM_THREADS_PER_BLOCK 256
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>







__global__ static void d_cross_correlation_obs(float* invalues, float* ans, int max_N) {
  const int idx = threadIdx.x + blockIdx.x*NUM_THREADS_PER_BLOCK;
  
  int cur_idx;
  int min_idx = max(0, idx - MAX_XCORR_OFFSET);
  int max_idx = min(max_N-1, idx + MAX_XCORR_OFFSET);

  ans[idx] = invalues[idx];
  
  for (cur_idx = min_idx ; cur_idx <= max_idx ; cur_idx++) {
    ans[idx] -= (invalues[cur_idx] / (MAX_XCORR_OFFSET * 2.0));
  }

  //ans[idx] = min_idx;
  //ans[idx] = idx;
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
void cross_correlation_obs(float* h_values, float* h_ans, int n) {
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

  d_cross_correlation_obs<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_values, d_ans, n);

  //check for errors.
  error = cudaGetLastError();

  if (error != cudaSuccess)
    printf ("Kernel error:%i\n",error);
  
  //copy answer back.
  cudaMemcpy(h_ans, d_ans, sizeof(float)*n, cudaMemcpyDeviceToHost);

}
 






#endif
