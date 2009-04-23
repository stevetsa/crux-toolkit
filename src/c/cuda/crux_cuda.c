#include "crux_cuda.h"

#include "scorer.h"

#include "cublas.h"

#include <math.h>

#include "mytimer.h"


BOOLEAN_T cudablas_initialized = FALSE;

#ifdef CRUX_USE_CUDA1
int current_size = -1;
float* d_cutemp = NULL;
float* d_observed = NULL;
float* d_theoretical = NULL;
#endif

#ifdef CRUX_USE_CUDA2
int spectrum_size = -1;
int max_spectrum_size = -1;

float* h_theoretical_matrix = NULL;


float* d_pThe_Matrix = NULL;
float* d_pObs = NULL;
float* d_pXCorrs = NULL;
#endif


struct my_timer timer2;
struct my_timer timer3;

BOOLEAN_T initialize_cudablas() {
  my_timer_reset(&timer2);
  my_timer_reset(&timer3);
  cublasStatus status;

  carp(CARP_FATAL,"initialize_cudablas(): begin");

  if (cudablas_initialized)
    return TRUE;
  else
    {
      //Initialize cudablas.
      status = cublasInit();

      if (status != CUBLAS_STATUS_SUCCESS) {
	carp(CARP_FATAL, "!!!! CUBLAS initialization error:%d\n", cublasGetError());
	return FALSE;
      }
      else {
	cudablas_initialized = TRUE;
	return TRUE;
      }
    }
}

BOOLEAN_T shutdown_cudablas() {
  cublasStatus status;

  printf("setVector\n");
  my_timer_print(&timer2);
  printf("dot product\n");
  my_timer_print(&timer3);

#ifdef CRUX_USE_CUDA1
  if (d_cutemp != NULL) cublasFree(d_cutemp);
  if (d_observed != NULL) cublasFree(d_observed);
  if (d_theoretical != NULL) cublasFree(d_theoretical);
#endif
#ifdef CRUX_USE_CUDA2
  if (d_pThe_Matrix != NULL)
    cublasFree(d_pThe_Matrix);
  if (d_pObs != NULL)
    cublasFree(d_pObs);
  if (d_pXCorrs != NULL)
    cublasFree(d_pXCorrs);
#endif
  if (!cudablas_initialized)
    return TRUE;
  else {
    //shutdown cudablas.
    status = cublasShutdown();
    if (status != CUBLAS_STATUS_SUCCESS) {
      carp(CARP_FATAL, "!!!! CUBLAS shutdown error:%d\n", cublasGetError());
      return FALSE;
    }
    else {
      cudablas_initialized = FALSE;
      return TRUE;
    }
  }
}


BOOLEAN_T is_cudablas_initialized() {
  return cudablas_initialized; 
}

#ifdef CRUX_USE_CUDA1


BOOLEAN_T checkSpace(int size) {

  //so if the current allocated space is enough for the vector,
  //then we have enough.
  if (size <= current_size)
    return TRUE;
  else {
    //carp(CARP_FATAL," allocating %d space",size);
    //free the previous vectors
    if (d_theoretical != NULL) cublasFree(d_theoretical);
    CUBLASEXEC(cublasAlloc(size, sizeof(float), (void**)&d_theoretical),"alloc(theo)");

    if (d_observed != NULL) cublasFree(d_observed);
    CUBLASEXEC(cublasAlloc(size, sizeof(float), (void**)&d_observed),"alloc(obs)");

    if (d_cutemp != NULL) cublasFree(d_cutemp);
    CUBLASEXEC(cublasAlloc(size, sizeof(float), (void**)&d_cutemp),"alloc(temp)");

    current_size = size;
    return TRUE;
  }
}


void cuda_set_observed(float* raw_values, int n, int num_regions, 
		       int region_selector, int max_offset) {
  //my_timer_start(&timer2);
  //check space.
  checkSpace(n);
  //push raw data into vector.
  cublasSetVector(n, sizeof(float), raw_values, 1, d_cutemp, 1);
  //run the kernel.
  d_cuda_sqrt_max_normalize_and_cc2(d_cutemp, d_observed, n, num_regions, region_selector, max_offset);
  //my_timer_stop(&timer2);
}


/**
 * Uses an iterative cross correlation using CUDA Blas dot product fxn.
 *
 *\return the final cross correlation score between the observed and the
 *theoretical spectra
 */
float cross_correlation_cuda(
  float* h_theoretical, ///< the theoretical spectrum to score against the observed spectrum -in
  int size
  ) {
  float score_at_zero = 0;

  //copy theoretical to the device.
  my_timer_start(&timer2);
  score_at_zero = *h_theoretical;
  //CUBLASEXEC(cublasSetVector(size, sizeof(float), h_theoretical, 1, d_theoretical, 1),
  //	     "cross_correlation_cuda: setVector(theoretical)");
  my_timer_stop(&timer2);
  //observed should have already been set by cuda_set_observed.

  //get dot product.

  my_timer_start(&timer3);
  score_at_zero = cublasSdot(size, d_theoretical, 1, d_observed, 1); 

  my_timer_stop(&timer3);

  return score_at_zero / 10000.0;
}
#endif /*CRUX_USE_CUDA1*/

#ifdef CRUX_USE_CUDA2


void cuda_realloc_float(int n,  float** ptr) {
  if (*ptr !=NULL)
    CUBLASEXEC(cublasFree(*ptr));
  CUBLASEXEC(cublasAlloc(n, sizeof(float), (void**)ptr));
}

void cuda_set_spectrum_size(int size) {
  //carp(CARP_FATAL,"cuda_set_spectrum_size: start %d",size);
  if (size != spectrum_size) {
    carp(CARP_FATAL,"Updating size %d to %d",spectrum_size, size);
    if (max_spectrum_size < size) {
      int n = size;
      int m = cuda_get_max_theoreticalN(n);
      int n2 = n * m;
      
      carp(CARP_FATAL,"Allocating space n:%i n2:%d",n, n2);
      //allocate space for the theoretical matrix, the observed vector, and the result vector.
      cuda_realloc_float(n2, &d_pThe_Matrix);
      cuda_realloc_float(n, &d_pObs);
      cuda_realloc_float(m, &d_pXCorrs);
      max_spectrum_size = size;
      } 
    spectrum_size = size;
   }
  //carp(CARP_FATAL,"cuda_set_spectrum_size: done");
}

int cuda_get_max_theoretical(void) {
  return cuda_get_max_theoreticalN(spectrum_size);
}
int cuda_get_max_theoreticalN(int ssize) {
  //carp(CARP_FATAL,"ssize:%d",ssize);
  if (ssize < 1) return 0;

  //TODO: calculate based upon card memory.
  int size = 1 * 1024 * 1024; //1 MB of memory.
  int fsize = size / sizeof(float); //number of floats.

  //m*n + n + m = S
  //m - number of theoretical.
  //n - spectrum size.
  //m*n - matrix of theoreticals.
  //n - storage space for observed spectrum.
  //m - storage space for calculated xcorrs.

  //carp(CARP_FATAL,"fsize:%d ssize:%d sizeof(float):%d",fsize,ssize,sizeof(float));
  

  int ans = (fsize - ssize) / (ssize + 1);
  
  //carp(CARP_FATAL,"number of theoretical:%d",ans);

  return ans;
}


void cuda_set_theoretical(float* h_pThe, int index) {
  float* h_ptr = h_theoretical_matrix + (spectrum_size * index);  
  memcpy(h_ptr,h_pThe,sizeof(float)*spectrum_size);
}

/*
void cuda_set_observed(float* h_pObs) {
  cublasStatus status;
  status = cublasSetVector(spectrum_size, sizeof(float), h_pObs, 1, d_pObs, 1);
  if (status != CUBLAS_STATUS_SUCCESS) {
    carp(CARP_FATAL,"!!!! Error setting obsered. code:%i",status);
  }
}
*/

void cuda_calculate_xcorrs(float* xcorrs) {
  cuda_calculate_xcorrsN(xcorrs, cuda_get_max_theoretical());
}

void cuda_calculate_xcorrsN(float* xcorrs, int nthe) {
  cublasStatus status;

  //copy the theoretical to the device.
  status = cublasSetVector(d_theoretical_matrix, h_theoretical_matrix, 
			   nthe*spectrum_size*sizeof(float), cudaMemcpyHostToDevice);


  /*
    trans specifies op(A). If trans == 'N' or 'n', .
    If trans == 'T', 't', 'C', or 'c', .
    m specifies the number of rows of matrix A; m must be at least zero.
    n specifies the number of columns of matrix A; n must be at least zero.
    alpha single-precision scalar multiplier applied to op(A).
    A single-precision array of dimensions (lda, n) if trans == 'N' or
    'n', of dimensions (lda, m) otherwise; lda must be at least
    max(1, m) if trans == 'N' or 'n' and at least max(1, n) otherwise.
    lda leading dimension of two-dimensional array used to store matrix A.
    x single-precision array of length at least if
    trans == 'N' or 'n', else at least .
    incx specifies the storage spacing for elements of x; incx must not be zero.
    beta single-precision scalar multiplier applied to vector y. If beta is zero, y
    is not read.
  */
  char trans = 'T'; //transpose to get dot products.
  int m = spectrum_size; //rows is the spectrum_size.
  int n = nthe; //columns are the number of theoretical spectra.
  float alpha = 1.0 / 10000.0; //for the wierd division part of the cross-correlation calculation.
  float* A = d_theoretical_matrix;
  int lda = m; //this has to equal the number of rows???
  float* x = d_observed;
  int incx = 1;
  float beta = 0;
  float* y = d_xcorrs;
  int incy = 1;

  //carp(CARP_FATAL, "CRUX_USE_CUDA2:Multiplying");
  //do the calculation.
  cublasSgemv (trans, m, n, alpha, A, lda, x,
	       incx, beta, y, incy);
  //check the error.
   status = cublasGetError();
   if (status != CUBLAS_STATUS_SUCCESS) {
     carp(CARP_FATAL,"!!!! kernel execution error. %i",status);
    }
    
   //carp(CARP_FATAL, "CRUX_USE_CUDA2:recovering result");
  /* Read the result back */
  status = cublasGetVector(nthe, sizeof(float), d_pXCorrs, 1, xcorrs, 1); 

}

#endif /*CRUX_USE_CUDA2*/

#ifdef CRUX_USE_CUDA3

typedef struct theoretical_spectra_t {
   int charge;
   int spectrum_size;
   float min_mass;
   float max_mass;
   int max_size;
   int current_size;
   float* mass;
  int* count;
   float* d_matrix;
}THEORETICAL_SPECTRA_T;




int calculateSpectrumSizeMass(float m1) {
  float experimental_mass_cut_off = m1 + 50 + 3 /*for +/-mass search*/;
  int ans = 512;
  // set max_mz and malloc space for the observed intensity array
  if(experimental_mass_cut_off > 512){
    int x = (int)experimental_mass_cut_off / 1024;
    float y = experimental_mass_cut_off - (1024 * x);
    ans = x * 1024;

    if(y > 0){
      ans += 1024;
    }
  }  
  return ans;
}

int calculateSpectrumSizeMZ(float mz1, int charge) {
  return calculateSpectrumSizeMass(mz1 * (float)charge);
}


#define max_ssize 8
THEORETICAL_SPECTRA_T theoretical_data[3][max_ssize]; /*charge,spectrum_size*/

int spectrum_sizes[] = {512, 1024, 2048, 3072, 4096, 5120, 6144, 7168};




void crux_cuda_init() {
  int s_index;
  int c_index;
  for (s_index = 0;s_index < max_ssize;s_index++)
    for (c_index = 0;c_index < 3;c_index++) {
      theoretical_data[c_index][s_index].charge = c_index + 1;
      theoretical_data[c_index][s_index].spectrum_size = spectrum_sizes[s_index];
      theoretical_data[c_index][s_index].min_mass = 10000;
      theoretical_data[c_index][s_index].max_mass = -1;
      theoretical_data[c_index][s_index].max_size = 80000;
      theoretical_data[c_index][s_index].current_size = 0;
      theoretical_data[c_index][s_index].mass = malloc(sizeof(float)*
						      theoretical_data[c_index][s_index].max_size);
      theoretical_data[c_index][s_index].count = malloc(sizeof(int)*
							theoretical_data[c_index][s_index].max_size);

      //TODO Allocate space on the device.
      //theoretical_data[s_index][c_index].d_matrix
    }
  printStats();
}

float min(float a, float b) {
  return a < b ? a : b;
}
float max(float a, float b) {
  return a > b ? a : b;
}

void generateMatrices(MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator,int charge) {
 
  THEORETICAL_SPECTRA_T *data;

  while( modified_peptides_iterator_has_next(peptide_iterator)){
    // get peptide
    PEPTIDE_T* peptide = modified_peptides_iterator_next(peptide_iterator);
    
    //assume the peptides are sorted by mass.
    
    //first we need to spectrum_size and the spectrum_index.
    //retrieve the mass.
    float mass = get_peptide_peptide_mass(peptide);
    int ssize = calculateSpectrumSizeMass(mass);
    int sindex = getSIndex(ssize);



    //if we are handling this spectrum size and there is still space
    //on the device, create a theoretical
    if (sindex != -1) {
      
      data = &theoretical_data[charge-1][sindex];
      //okay, we are going to create a theoretical spectrum for this peptide.

      //carp(CARP_FATAL,"updating");
      
      float mz = mass / (float)charge;
      if (mz >= 400 && mz <=1400) {
	data -> mass[data -> current_size] = mass;
	data -> count[data -> current_size] = 0;
	data -> max_mass = max(data -> max_mass, mass);
	data -> min_mass = min(data -> min_mass, mass);
	data -> current_size++;
	/*
	  carp(CARP_FATAL, "c:%i s:%i mz1:%f max:%f min:%f size:%i",
	  charge,
	  sindex,
	  mass,
	  data -> max_mass,
	  data -> min_mass,
	  data -> current_size);
	*/
      }
    }
    else {
      carp(CARP_FATAL,"spectrum size not handled:%i",ssize);
    }
  }

  

  carp(CARP_FATAL,"done!");
}


void printStats() {
  int sindex;
  int cindex;
  int mindex;
  THEORETICAL_SPECTRA_T *ptr;
  for (cindex=0;cindex<3;cindex++) {
    for (sindex=0;sindex<max_ssize;sindex++) {
      ptr = &theoretical_data[cindex][sindex];
      carp(CARP_FATAL,"threoretical[%i][%i]",cindex,sindex);
      carp(CARP_FATAL,"charge:%i",ptr -> charge);
      carp(CARP_FATAL,"spectrum_size:%i",ptr -> spectrum_size);
      carp(CARP_FATAL,"min_mass:%f", ptr -> min_mass);
      carp(CARP_FATAL,"max_mass:%f", ptr -> max_mass);
      carp(CARP_FATAL,"max_size:%i", ptr -> max_size);
      carp(CARP_FATAL,"current:%i",ptr -> current_size);
      for(mindex=0;mindex<ptr -> current_size;mindex++) {
	if (ptr -> count[mindex] != 0)
	  carp(CARP_FATAL,"mass:%f count:%i",
	       ptr -> mass[mindex],
	       ptr -> count[mindex]);
      }
      carp(CARP_FATAL,"===================");
    }
  }
}


void writeStats() {
  int sindex;
  int cindex;
  int mindex;

  FILE*fp;
  char buf[50];

  THEORETICAL_SPECTRA_T *ptr;
  for (cindex=0;cindex<3;cindex++) {
    for (sindex=0;sindex<max_ssize;sindex++) {
      ptr = &theoretical_data[cindex][sindex];
      if (ptr -> current_size != 0) {
	sprintf(buf,"theo_%i_%i",ptr -> charge, ptr -> spectrum_size);
	fp = fopen(buf, "w");
	for(mindex=0;mindex<ptr -> current_size;mindex++) {
	  fprintf(fp,"%f\t%i\n",ptr -> mass[mindex],ptr -> count[mindex]);
	}
	fclose(fp);
      }
    }
  }
}


void updateStatCount(float mass, int charge) {
  int ssize = calculateSpectrumSizeMass(mass);
  int sindex = getSIndex(ssize);
  int cindex = charge-1;
  int i;
  THEORETICAL_SPECTRA_T* ptr = &theoretical_data[cindex][sindex];
  
  double min_mass = mass - 3.0;
  double max_mass = mass + 3.0;

  if (min_mass < ptr -> mass[0])
    carp(CARP_FATAL,"min:%f min:%f",min_mass,ptr -> mass[0]);

  for (i=0;i<ptr -> current_size;i++) {
    if (ptr -> mass[i] < min_mass) continue;
    else if (ptr -> mass[i] > max_mass) break;
    else {
      ptr -> count[i]++;
    }
  }


}



int getSIndex(int ssize) {
  int ans = -1;
  int i;
  for (i=0;i<max_ssize;i++) {
    if (spectrum_sizes[i] == ssize) {
      ans = i;
      break;
    }
  }
  return ans;
}

#endif /*CRUX_USE_CUDA3*/
