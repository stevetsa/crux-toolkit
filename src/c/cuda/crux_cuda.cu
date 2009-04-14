
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

__global__ static void reduction_odd(float* a, int stride, int N) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int n = stride*2;
  int index = i*n;
  if (index < N) {
    if (index+stride != N) {
      int index = i * n;
      //leave the last element alone, it is the max in the first stride.
      a[index] = max(a[index], a[index+stride]);
    }
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
  cudaError error;
  float* d_values;
  float* d_ans;
  int size_n = sizeof(float) * n;


  error = cudaMalloc((void**)&d_values, size_n);
  error = cudaMalloc((void**)&d_ans, sizeof(float));
  
  error = cudaMemcpy(d_values, h_values, size_n,cudaMemcpyHostToDevice);
  
  calcMax(d_values, n, d_ans);

  error = cudaMemcpy(&ans, d_ans, sizeof(float), cudaMemcpyDeviceToHost);

  cudaFree(d_values);
  cudaFree(d_ans);
  
  return ans;

}

int isPowerOfTwo(int n)
{
    return (n) && !(n & (n - 1)); //this checks if the integer n is a power of two or not
}

void calcMaxP2(float* d_values, int n, float* d_ans);
void calcMaxNP2(float* d_values, int n, float* d_ans);

void calcMax(float* d_values, int n, float*d_ans) {
  


  if (isPowerOfTwo(n)) {
    calcMaxP2(d_values,n,d_ans);
  }
  else {
    calcMaxNP2(d_values,n,d_ans);
  }
}

void calcMaxP2(float* d_values, int n,float* d_ans) {
  int stride;
  cudaError error;
  int i;
  int size_n = n * sizeof(float);
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

  printf("n:%d num_threads:%d num_blocks:%d\n",n,num_threads, num_blocks);

  for(stride=1;stride<=n/2;stride*=2) {
    reduction<<<num_blocks, num_threads>>>(d_values, stride, n);
    cudaThreadSynchronize();
    error = cudaGetLastError();
    if (error != cudaSuccess)
      printf("reduction error %s stride:%i n:%i\n",cudaGetErrorString(error),stride, n);
    

  }
  cudaMemcpy(d_ans, 
	     d_values, 
	     sizeof(float), 
	     cudaMemcpyDeviceToDevice);
}


void calcMaxNP2(float* d_values, int n, float* d_ans){
  printf("Inside calcMaxNP2\n");
  int n2 = (int)pow(2,ceil(log2((float)n)));
  float* d_temp;
  cudaError error;
  
  printf("n:%i n2:%i\n",n, n2);

  cudaMalloc((void**)&d_temp, n2*sizeof(float));
  cudaMemcpy(d_temp,d_values,n*sizeof(float),cudaMemcpyDeviceToDevice);
  cudaMemset(d_temp+n,0,n2-n);
  calcMaxP2(d_temp, n2, d_ans);
  cudaFree(d_temp);

}

void calcMax2(float* d_values, int start, int end, float* d_ans) {
  calcMax(d_values+start,end-start,d_ans);
}


void cuda_sqrt_max_normalize_and_cc(float* h_values, int n, int num_regions, int region_selector, int max_offset) {
  float *d_values;
  float *d_max_per_region;
  float *d_ans;
  float *h_max_per_region;
  int i;
  int region;

  cudaError error;
  
  int size_n = n * sizeof(float);
  int size_reg = num_regions * sizeof(float);

  h_max_per_region = new float[num_regions];
  if (h_max_per_region == NULL) {
    printf("Error allocating h_max_per_region\n");
    exit(-1);
  }
  memset(h_max_per_region, 0, size_reg);

  //printf("allocating memory %i %i\n",size_n, size_reg);

  error = cudaMalloc((void**)&d_values, size_n);
  if (error != cudaSuccess)
    printf("allocating d_values memory error:%i\n",error);
  error = cudaMalloc((void**)&d_max_per_region, size_reg);
  if (error != cudaSuccess)
    printf("allocating d_max_per_region error:%i\n",error);
  error = cudaMalloc((void**)&d_ans, size_n); 
  if (error != cudaSuccess)
    printf("allocating d_ans memory error:%i\n",error);

  //printf("Copying data\n");
  error = cudaMemcpy(d_values, h_values, size_n, cudaMemcpyHostToDevice);
  if (error != cudaSuccess)
    printf("error copying h_values => d_Values:%i\n",error);

  //error = cudaMemcpy(d_max_per_region, h_max_per_region, size_reg, cudaMemcpyHostToDevice); 
 
  int num_blocks = n / NUM_THREADS_PER_BLOCK; 

  //printf("Executing d_sqrt_and_max_region\n");
  d_sqrt<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_values);
 
  error = cudaGetLastError(); 
 
  if (error != cudaSuccess) 
    printf("d_sqrt: error %i\n",error); 
 
  error = cudaMemcpy(d_ans, d_values, size_n, cudaMemcpyDeviceToDevice);

  for(region=0;region<num_regions;region++) {
    int start = region_selector * region;
    int end = min(n, start + region_selector);
    printf("calculating max for region %i %i %i\n",region, start, end);
    calcMax2(d_ans, start, end, d_max_per_region+region);
  }


  //printf("Copying max_regions\n");

  cudaMemcpy(h_max_per_region,d_max_per_region, size_reg, cudaMemcpyDeviceToHost);

  for (i=0;i<num_regions;i++)
    printf("cuda max[%i]:%f\n",i,h_max_per_region[i]);

  //printf("Normalizing regions\n");

  d_normalize_each_region<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_values, d_max_per_region,  
                                                                 n, num_regions, region_selector); 
  error = cudaGetLastError(); 
  
  if (error != cudaSuccess)
    printf("d_normalize_each_region: error %i\n",error);

  d_cross_correlation_obs<<<num_blocks, NUM_THREADS_PER_BLOCK>>>(d_values, d_ans, n, max_offset); 
  error = cudaGetLastError(); 

  if (error != cudaSuccess)
    printf("d_cross_correlation_obs: error %i\n", error);
 
  error = cudaMemcpy(h_values, d_ans, size_n, cudaMemcpyDeviceToHost); 
  
  //printf("cuda do cleanup\n");
  free(h_max_per_region);
  cudaFree(d_values);
  cudaFree(d_max_per_region);
  cudaFree(d_ans);
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
