#include "scorer.h"
#include "utils.h"
#ifdef CRUX_USE_CUDA3
#include "modified_peptides_iterator.h"
#endif


BOOLEAN_T initialize_cudablas();

BOOLEAN_T shutdown_cudablas();
BOOLEAN_T is_cudablas_initialized();

#ifdef CRUX_USE_CUDA
float cross_correlation_cuda(
  float* theoretical, ///< the theoretical spectrum to score against the observed spectrum -in
  float* observed,
  int size
  );
#endif


#ifdef CRUX_USE_CUDA2
//Sets the size of the vector for calculating the xcorrs.
void cuda_set_spectrum_size(int size);

//returns the maximum number of theoretical spectra that can be calculated
//by the card for a set spectrum_size.
int cuda_get_max_theoretical();

//sets the intensity array for a theoretical spectrum.
void cuda_set_theoretical(float* h_pThe, int index);

//sets the intensity array for the observed spectrum.
void cuda_set_observed(float*h_pObs);

//runs the calculation and copies the results to xcorrs.
void cuda_calculate_xcorrs(float* xcorrs);
void cuda_calculate_xcorrsN(float* xcorrs, int nthe);
#endif

#ifdef CRUX_USE_CUDA3
void crux_cuda_init();
void generateMatrices(MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator,int charge);
int getSIndex(int ssize);
void printStats();
void writeStats();
void updateStatCount(float mass, int charge);
#endif
