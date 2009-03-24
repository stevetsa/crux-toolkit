#include "scorer.h"
#include "utils.h"

BOOLEAN_T initialize_cudablas();

BOOLEAN_T shutdown_cudablas();
BOOLEAN_T is_cudablas_initialized();

float cross_correlation_cuda(
  float* theoretical, ///< the theoretical spectrum to score against the observed spectrum -in
  float* observed,
  int size
  );

