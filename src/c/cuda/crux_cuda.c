#include "crux_cuda.h"

#include "scorer.h"

#include "cublas.h"

BOOLEAN_T cudablas_initialized = FALSE;

float *d_pThe = NULL;
float *d_pObs = NULL;
int current_size = -1;
float *h_pObs_current = NULL;


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

  //check to see if vectors are allocated.
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
