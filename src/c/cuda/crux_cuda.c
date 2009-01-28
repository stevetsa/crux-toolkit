#include "crux_cuda.h"

#include "scorer.h"

#include "cublas.h"


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
  )
{
  float score_at_zero = 0;
  //int i;
  //pointers to device memory
  float *d_pThe;
  float *d_pObs;

  //Initialize cudablas.
  cublasStatus status = cublasInit();

  if (status != CUBLAS_STATUS_SUCCESS) 
    {
      carp(CARP_FATAL, "!!!! CUBLAS initialization error\n");
      return 0.0;
    }



  //allocate space for both vectors.
  status = cublasAlloc(size, sizeof(float), (void**)&d_pThe);
   if (status != CUBLAS_STATUS_SUCCESS) 
    {
      carp(CARP_FATAL, "!!!! CUBLAS alloc theoretical error\n");
      return 0.0;
    }

  status = cublasAlloc(size, sizeof(float), (void**)&d_pObs);
  if (status != CUBLAS_STATUS_SUCCESS) 
    {
      carp(CARP_FATAL, "!!!! CUBLAS alloc observed error\n");
      return 0.0;
    }

  /*
  for (i=0;i<size;i++) {
    if (observed[i] != 0) {
      carp(CARP_FATAL,"theoretical[%i]=%f",i, theoretical[i]);
      carp(CARP_FATAL,"observed[%i]=%f",i, observed[i]);
    }
  }
  */
  //copy vectors to device.
  status = cublasSetVector(size, sizeof(float), theoretical, 1, d_pThe, 1);
  if (status != CUBLAS_STATUS_SUCCESS) 
    {
      carp(CARP_FATAL, "!!!! CUBLAS set theoretical error\n");
      return 0.0;
    }
  status = cublasSetVector(size, sizeof(float), observed, 1, d_pObs, 1);
  if (status != CUBLAS_STATUS_SUCCESS) 
    {
      carp(CARP_FATAL, "!!!! CUBLAS set observed error\n");
      return 0.0;
    }
  

  //get dot product.
  score_at_zero = cublasSdot(size, d_pThe, 1, d_pObs, 1); 
  
  //free memory on device
  cublasFree(d_pThe);
  cublasFree(d_pObs);

  //shutdown cudablas.
  status = cublasShutdown();

  //carp(CARP_FATAL,"score_at_zero:%f",score_at_zero);

  return score_at_zero / 10000.0;
}
