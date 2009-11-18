#ifndef CRUX_SCORER_H_
#define CRUX_SCORER_H_

extern "C" {
#include "scorer.h"

/**
 * create the intensity arrays for both observed and theoretical spectrum
 * SCORER must have been created for XCORR type
 * \returns TRUE if successful, else FLASE
 */
BOOLEAN_T create_intensity_array_xcorr(
  SPECTRUM_T* spectrum,    ///< the spectrum to score(observed) -in
  SCORER_T* scorer,        ///< the scorer object -in/out
  int charge               ///< the peptide charge -in 
  );

/**
 * Uses an iterative cross correlation
 *
 *\return the final cross correlation score between the observed and the
 *theoretical spectra
 */
FLOAT_T cross_correlation(
  SCORER_T* scorer,  ///< the scorer object that contains observed spectrum -in
  FLOAT_T* theoretical ///< the theoretical spectrum to score against the observed spectrum -in
  );


/**
 * adds the intensity at add_idx
 * if, there already exist a peak at the index, only overwrite if
 * intensity is larger than the existing peak.
 */
void add_intensity(
  FLOAT_T* intensity_array, ///< the intensity array to add intensity at index add_idx -out
  int add_idx,            ///< the idex to add the intensity -in
  FLOAT_T intensity         ///< the intensity to add -in
  );

}

#endif
