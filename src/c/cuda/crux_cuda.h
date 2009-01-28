#include "scorer.h"

float cross_correlation_cuda(
  float* theoretical, ///< the theoretical spectrum to score against the observed spectrum -in
  float* observed,
  int size
  );

