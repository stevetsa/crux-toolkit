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

  void d_cuda_sqrt_max_normalize_and_cc(float* d_in, int n, int num_regions, 
					int region_selector, int max_offset);

  void d_cuda_sqrt_max_normalize_and_cc2(float* d_in, float* d_out, int n, int num_regions,
					 int region_selector, int max_offset);

  void cuda_sqrt_max_normalize_and_cc(float* h_values, int n, int num_regions, 
				      int region_selector, int max_offset);
 
  
  void cuda_sort(float* d_values, short* d_index, int n);

  void cuda_float_malloc(float** a, int current_size, int size);
  
  void crux_cuda_initialize();
  void crux_cuda_shutdown();

#ifdef CUDA_NVCC
}
#endif
/*
#define CUDAEXEC(x,s) {cudaError _status; \
                      _status = x; \
                      if (_status != cudaSuccess) \
                      {printf("cuda error: %s: %s\n", cudaGetErrorString(_status),s);}}
*/
#define CUDAEXEC(x, s) x;
