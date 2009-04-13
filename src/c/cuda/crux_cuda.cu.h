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
  
  void cuda_normalize_each_region(float* h_values, float* max_per_region, int n, int num_regions, int region_selector);
  void cross_correlation_obs(float* h_values, float*h_ans, int n, int max_offset);
  


#ifdef CUDA_NVCC
}
#endif

