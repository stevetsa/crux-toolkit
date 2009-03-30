#include "crux_cuda.h"

#include "scorer.h"

#include "cublas.h"

BOOLEAN_T cudablas_initialized = FALSE;

#ifdef CRUX_USE_CUDA2
int spectrum_size = -1;
float* d_pThe_Matrix = NULL;
float* d_pObs = NULL;
float* d_pXCorrs = NULL;
#endif


BOOLEAN_T initialize_cudablas() {

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
#ifdef CRUX_USE_CUDA1
  if (d_pThe != NULL) {
      //free memory on device
      cublasFree(d_pThe);
      d_pThe = NULL;
    }
  if (d_pObs != NULL) {
    //free memory on device
    cublasFree(d_pObs);
    d_pObs = NULL;
  }
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
float *d_pThe = NULL;
float *d_pObs = NULL;
int current_size = -1;
float *h_pObs_current = NULL;

BOOLEAN_T checkSpace(int size) {
  cublasStatus status;

  //so if the current allocated space is enough for the vector,
  //then we have enough.
  if (size <= current_size)
    return TRUE;
  else {
    //carp(CARP_FATAL," allocating %d space",size);
    //free the previous vectors
    if (d_pThe != NULL)
      cublasFree(d_pThe);
    if (d_pObs != NULL)
      cublasFree(d_pObs);
    
    //allocate space for both vectors.
    status = cublasAlloc(size, sizeof(float), (void**)&d_pThe);
    if (status != CUBLAS_STATUS_SUCCESS) 
      {
	carp(CARP_FATAL, "!!!! CUBLAS alloc theoretical error\n");
	return FALSE;
      }
    
    status = cublasAlloc(size, sizeof(float), (void**)&d_pObs);
    if (status != CUBLAS_STATUS_SUCCESS) 
      {
	carp(CARP_FATAL, "!!!! CUBLAS alloc observed error\n");
	return FALSE;
      }
    //update the size.
    current_size = size;
    return TRUE;
  }
}


/**
 * Uses an iterative cross correlation using CUDA Blas dot product fxn.
 *
 *\return the final cross correlation score between the observed and the
 *theoretical spectra
 */
float cross_correlation_cuda(
  float* theoretical, ///< the theoretical spectrum to score against the observed spectrum -in
  float* observed,
  int size
  ) {
  float score_at_zero = 0;
  
  cublasStatus status;
  
  BOOLEAN_T size_changed = size != current_size;

  //check to see if there is enough allocated space on the device.
  checkSpace(size);

  //copy theoretical to the device.
  status = cublasSetVector(size, sizeof(float), theoretical, 1, d_pThe, 1);
  if (status != CUBLAS_STATUS_SUCCESS) 
    {
      carp(CARP_FATAL, "!!!! CUBLAS set theoretical error\n");
      return 0.0;
    }

  //assume that the observed won't change for a great number of
  //theoretical spectra, i.e. only copy the observed to the 
  //device if the observed changes.
  if (h_pObs_current != observed || size_changed) {
    status = cublasSetVector(size, sizeof(float), observed, 1, d_pObs, 1);
    if (status != CUBLAS_STATUS_SUCCESS) 
      {
	carp(CARP_FATAL, "!!!! CUBLAS set observed error\n");
	return 0.0;
      }
    h_pObs_current = observed;
  }
  else {
    //carp(CARP_FATAL,"not copying observed");
  }
  //get dot product.
  score_at_zero = cublasSdot(size, d_pThe, 1, d_pObs, 1); 

  return score_at_zero / 10000.0;
}
#endif /*CRUX_USE_CUDA*/

//float *h_temp = NULL;
//float *h_temp2 = NULL;
#ifdef CRUX_USE_CUDA2
void cuda_set_spectrum_size(int size) {
  cublasStatus status;

  if (size != spectrum_size) {
    //carp(CARP_FATAL, "cuda_set_spectrum_size():%d",size);
    spectrum_size = size;
    int n = spectrum_size;
    int n2 = spectrum_size * cuda_get_max_theoretical();
    //allocate space for the theoretical matrix, the observed vector, and the result vector.
    if (d_pThe_Matrix != NULL) cublasFree(d_pThe_Matrix);
    status = cublasAlloc(n2, sizeof(float), (void**)&d_pThe_Matrix);

    if (d_pObs != NULL) cublasFree(d_pObs);
    status = cublasAlloc(n, sizeof(float), (void**)&d_pObs);
    /*
    if (h_temp != NULL) free(h_temp);
    h_temp = malloc(sizeof(float)*n);
    
    if (h_temp2 != NULL) free(h_temp2);
    h_temp2 = malloc(sizeof(float)*n);
    */

    if (d_pXCorrs != NULL) cublasFree(d_pXCorrs);
    status = cublasAlloc(cuda_get_max_theoretical(), sizeof(float), (void**)&d_pXCorrs);
    if (status != CUBLAS_STATUS_SUCCESS) {
      carp(CARP_FATAL,"!!!! Error allocating d_pXCorrs. %i",status);
    }

  }

}

int cuda_get_max_theoretical() {
  //TODO: calculate based upon card memory.
  return 10;
}


void cuda_set_theoretical(float* h_pThe, int index) {
  cublasStatus status;
  /*
  int i;
  for (i=0;i<spectrum_size;i++)
    carp(CARP_FATAL, "Cuda_set_theoretical[%i]=%f",i,h_pThe[i]);
  */
  //carp(CARP_FATAL, "cuda_set_theoretical:index:%d",(index+1));
  status = cublasSetVector(spectrum_size, sizeof(float), h_pThe, 1, d_pThe_Matrix, index+1);
}

void cuda_set_observed(float* h_pObs) {
  cublasStatus status;
  /*
  int i;
  for (i=0;i<spectrum_size;i++)
    carp(CARP_FATAL, "cuda_set_observed[%i]=%f",i,h_pObs[i]);
  */

  status = cublasSetVector(spectrum_size, sizeof(float), h_pObs, 1, d_pObs, 1);
}

void cuda_calculate_xcorrs(float* xcorrs) {
  cuda_calculate_xcorrsN(xcorrs, cuda_get_max_theoretical());
}

void cuda_calculate_xcorrsN(float* xcorrs, int nthe) {
  cublasStatus status;
  //carp(CARP_FATAL, "cuda_calculate_xcorrsN:%d",nthe);

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
  float* A = d_pThe_Matrix;
  int lda = m; //this has to equal the number of rows???
  float* x = d_pObs;
  int incx = 1;
  float beta = 0;
  float* y = d_pXCorrs;
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
