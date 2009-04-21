
#ifndef _CRUX_CUDA_CU_
#define _CRUX_CUDA_CU_

#include "crux_cuda.cu.h"

#include <stdio.h>

#define NUM_THREADS_PER_BLOCK 256

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
//#include "carp.h"


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

void calcMax(float* d_values, int n, float*d_ans);
float d_max(float* h_values, int n) {

  float ans;
  float* d_values;
  float* d_ans;
  int size_n = sizeof(float) * n;


  CUDAEXEC(cudaMalloc((void**)&d_values, size_n),"alloc d_values");
  CUDAEXEC(cudaMalloc((void**)&d_ans, sizeof(float)),"alloc d_ans");
  CUDAEXEC(cudaMemcpy(d_values, h_values, size_n,cudaMemcpyHostToDevice),"d_val -> h_val");
  calcMax(d_values, n, d_ans);

  CUDAEXEC(cudaMemcpy(&ans, d_ans, sizeof(float), cudaMemcpyDeviceToHost),"d_ans -> ans");

  CUDAEXEC(cudaFree(d_values),"Free d_values");
  CUDAEXEC(cudaFree(d_ans),"Free d_ans");
  
  return ans;

}

int isPowerOfTwo(int n)
{
    return (n) && !(n & (n - 1)); //this checks if the integer n is a power of two or not
}

void calcMaxP2(float* d_values, int n, float* d_ans);
void calcMaxNP2(float* d_values, int n, float* d_ans);
/****************************************************
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
void calcMax(float* d_values, int n, float*d_ans) {
  if (isPowerOfTwo(n)) {
    calcMaxP2(d_values,n,d_ans);
  }
  else {
    calcMaxNP2(d_values,n,d_ans);
  }
}


/******************************************************
 *calcMaxP2: calculates the max in an array that has
 *a power of 2 elements
 ******************************************************/
void calcMaxP2(float* d_values, int n,float* d_ans) {
  int stride;
  cudaError error;
  int num_threads;
  int num_blocks;
  if (n < NUM_THREADS_PER_BLOCK) {
    num_threads = n;
    num_blocks = 1;
  }
  else {
    num_threads = NUM_THREADS_PER_BLOCK;
    num_blocks = n / num_threads;
    if (n % num_threads !=0) num_blocks++;
  }

  //printf("n:%d num_threads:%d num_blocks:%d\n",n,num_threads, num_blocks);

  for(stride=1;stride<=n/2;stride*=2) {
    reduction<<<num_blocks, num_threads>>>(d_values, stride, n);
    error = cudaGetLastError();
    if (error != cudaSuccess)
      printf("reduction error %s stride:%i n:%i\n",cudaGetErrorString(error),stride, n);
    cudaThreadSynchronize();
    

  }
  CUDAEXEC(cudaMemcpy(d_ans, 
		      d_values, 
		      sizeof(float), 
		      cudaMemcpyDeviceToDevice),
	   "calcMaxP2:d_values -> d_ans");
}

/*********************************************************
 *calcMaxNP2: calculates the max in an array that does not
 *have power of 2 elements. It achieves this by copying to
 *a temporary array that is a power of 2.
 *********************************************************/
void calcMaxNP2(float* d_values, int n, float* d_ans){
  //printf("Inside calcMaxNP2\n");
  int n2 = (int)pow(2,ceil(log2((float)n)));
  float* d_temp;
  
  //printf("n:%i n2:%i\n",n, n2);

  CUDAEXEC(cudaMalloc((void**)&d_temp, n2*sizeof(float)),"alloc d_temp");
  CUDAEXEC(cudaMemset(d_temp,0,n2*sizeof(float)),"memset(d_temp,0)");
  cudaError _error = cudaMemcpy(d_temp,d_values,n*sizeof(float),cudaMemcpyDeviceToDevice);
  if (_error != cudaSuccess) {
    printf("calcMaxNP2 memcpy d_values -> d_temp error:%s\n",cudaGetErrorString(_error));
    printf("n:%i n2:%i\n",n, n2);
  }
  calcMaxP2(d_temp, n2, d_ans);
  CUDAEXEC(cudaFree(d_temp),"cudaFree");

}

void calcMax2(float* d_values, int start, int end, float* d_ans) {
  if (end-start < 1) {
    cudaMemset(d_ans, 0, sizeof(float));
    //printf("calcMax2: end:%i start:%i\n",end, start);
  }
  else {
    calcMax(d_values+start,end-start,d_ans);
  }
}

void cuda_sqrt_max_normalize_and_cc(float* h_values, int n, int num_regions, int region_selector, int max_offset) {
  float *d_values;
  float *d_max_per_region;
  float *d_ans;
  int i;
  int region;

  cudaError error;
  
  int size_n = n * sizeof(float);
  int size_reg = num_regions * sizeof(float);

  //printf("allocating memory %i %i\n",size_n, size_reg);

  CUDAEXEC(cudaMalloc((void**)&d_values, size_n),"alloc d_values");
  CUDAEXEC(cudaMalloc((void**)&d_max_per_region, size_reg), "alloc d_max_per_region");
  CUDAEXEC(cudaMalloc((void**)&d_ans, size_n),"alloc d_ans"); 
  CUDAEXEC(cudaMemcpy(d_values, h_values, 
		      size_n, cudaMemcpyHostToDevice),
	   "h_values -> d_values");
 
  int num_blocks = n / NUM_THREADS_PER_BLOCK; 

  //printf("Executing d_sqrt_and_max_region\n");
  d_sqrt<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_values); 
  error = cudaGetLastError(); 
  if (error != cudaSuccess) 
    printf("d_sqrt: error %s\n",cudaGetErrorString(error)); 
 
  CUDAEXEC(cudaMemcpy(d_ans, d_values, size_n, cudaMemcpyDeviceToDevice),"d_values -> d_ans");

  for(region=0;region<num_regions;region++) {
    int start = region_selector * region;
    int end = min(n, start + region_selector);
    

    calcMax2(d_ans, start, end, d_max_per_region+region);
  }

  d_normalize_each_region<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_values, d_max_per_region,  
                                                                 n, num_regions, region_selector); 
  error = cudaGetLastError(); 
  
  if (error != cudaSuccess) {
    printf("d_normalize_each_region: error %s\n",cudaGetErrorString(error));
    printf("n: %i nr: %i s: %i nb: %i ntb: %i\n",n,num_regions,region_selector, num_blocks, NUM_THREADS_PER_BLOCK);
  }
  d_cross_correlation_obs<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_values, d_ans, n, max_offset); 
  error = cudaGetLastError(); 

  if (error != cudaSuccess)
    printf("d_cross_correlation_obs: error %s\n", cudaGetErrorString(error));
 
  CUDAEXEC(cudaMemcpy(h_values, d_ans, size_n, cudaMemcpyDeviceToHost),"d_ans -> h_values"); 
  
  //printf("cuda do cleanup\n");
  CUDAEXEC(cudaFree(d_values), "free(d_values)");
  CUDAEXEC(cudaFree(d_max_per_region), "free(d_max_per_region)");
  CUDAEXEC(cudaFree(d_ans), "free(d_ans)");
  //free(d_values);
  //free(d_max_per_region);
  //free(d_ans);
}


void cuda_normalize_and_cc(float*h_values, float* max_per_region,
			   int n, int num_regions,
			   int region_selector, int max_offset) {
  float* d_values;
  float* d_max_per_region;
  float* d_ans;

  cudaError error = cudaMalloc((void**)&d_values, sizeof(float)*n);
  if (error != cudaSuccess)
    {
      printf("Error allocating d_values:%i",error);
    }
  
  error = cudaMalloc((void**)&d_max_per_region, sizeof(float)* num_regions);
  error = cudaMalloc((void**)&d_ans, sizeof(float)*n);

  
  error = cudaMemcpy(d_values, h_values, sizeof(float)*n, cudaMemcpyHostToDevice);
  error = cudaMemcpy(d_max_per_region, max_per_region, sizeof(float)*n, cudaMemcpyHostToDevice);

  int num_blocks = n / NUM_THREADS_PER_BLOCK;

  d_normalize_each_region<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_values, d_max_per_region, 
								 n, num_regions, region_selector);
  error = cudaGetLastError();

  d_cross_correlation_obs<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_values, d_ans, n, max_offset);
  error = cudaGetLastError();

  error = cudaMemcpy(h_values, d_ans, sizeof(float)*n, cudaMemcpyDeviceToHost);
  

  //Free memory, (maybe some of these don't need to be freed everytime (TODO)).
  error = cudaFree(d_values);
  error = cudaFree(d_max_per_region);
  error = cudaFree(d_ans);
}



void cuda_normalize_each_region(float* h_values, float* max_per_region, 
				int n, int num_regions, int region_selector) {
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
    printf("Error allocating d_values:%s",cudaGetErrorString(error));

  error = cudaMalloc((void**)&d_ans, sizeof(float)*n);
  if (error != cudaSuccess) 
    printf("Error allocating d_ans:%s",cudaGetErrorString(error));

  
  error = cudaMemcpy(d_values, h_values, sizeof(float)*n, cudaMemcpyHostToDevice);
  if (error != cudaSuccess)
    printf("Error copying values:%s",cudaGetErrorString(error));

  int num_blocks = n / NUM_THREADS_PER_BLOCK;

  d_cross_correlation_obs<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_values, d_ans, n, max_offset);

  //check for errors.
  error = cudaGetLastError();

  if (error != cudaSuccess)
    printf ("Kernel error:%s\n",cudaGetErrorString(error));
  
  //copy answer back.
  cudaMemcpy(h_ans, d_ans, sizeof(float)*n, cudaMemcpyDeviceToHost);

}
 






#endif
