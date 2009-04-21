/*
 *CUDA KERNEL Functions
 */





//CUDA's nvcc compiler compiles code as c++
//which mangles the names.  Wrapping the code
//in extern "C" {} prevents this, however, C
//doesn't handle this well.  So the hack solution
//is to use a define which is passed in whenever
//NVCC is compiling the code that references this
//header, otherwise, the code will not be wrapped.
#ifdef CUDA_NVCC
extern "C" {
#endif

  

  void cuda_sqrt_max_normalize_and_cc(float* h_values, int n, int num_regions, int region_selector, int max_offset);
 
  
  void cuda_normalize_and_cc(float*h_values, float* max_per_region,
			     int n, int num_regions,
			     int region_selector, int max_offset);

  void cuda_normalize_each_region(float* h_values, float* max_per_region, 
				  int n, int num_regions, 
				  int region_selector);

  void cross_correlation_obs(float* h_values, float*h_ans, 
			     int n, int max_offset);
  
  float d_max(float* h_values, int n);


#ifdef CUDA_NVCC
}
#endif

#define CUDAEXEC(x,s) {cudaError _status; \
                      _status = x; \
                      if (_status != cudaSuccess) \
                      {printf("cuda error: %s: %s\n", cudaGetErrorString(_status),s);}}

