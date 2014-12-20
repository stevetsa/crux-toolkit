/**
 * \file Scorer.h 
 * \brief object to score spectrum vs. spectrum or spectrum vs. scorer
 */

/*
 * AUTHOR: Chris Park
 * CREATE DATE: 9 Oct 2006
 * $Revision: 1.22 $
 *****************************************************************************/
#ifndef SCORER_H 
#define SCORER_H

#include <stdio.h>
#ifndef _MSC_VER
#include <dirent.h>
#endif
#include <string>
#ifdef _MSC_VER
#include "windirent.h"
#endif
#include "objects.h"
#include "Spectrum.h"
#include "Peptide.h"
#include "Ion.h"


/**
 * Macro for converting floating point to integers.
 */

#define INTEGERIZE(VALUE,BIN_SIZE,BIN_OFFSET) \
  ((int)( ( ( VALUE / BIN_SIZE ) + 1.0 ) - BIN_OFFSET ) )

class Scorer {

 protected:
  SCORER_TYPE_T type_; ///< The type of scorer
  double sp_beta_; ///< used for Sp: the beta variable 
  double sp_max_mz_; ///< used for Sp: the max mz for the intensity array
  int sp_b_y_ion_matched_; ///< The most recent ion_collection number of the b, y ion matched while scoring for SP
  int sp_b_y_ion_possible_; ///< The most recent ion_collection number of the b, y ion possible while scoring for SP
  double sp_b_y_ion_fraction_matched_; ///< The ratio of matched and possible.

  double* intensity_array_; ///< used for Sp: the intensity array, which can be indexed using the m/z
  double max_intensity_; ///< the max intensity in the intensity array
  bool initialized_; ///< has the scorer been initialized?
  int last_idx_; ///< the last index in the array, the data size of the array

  double bin_width_; ///< width of the bins to use for arrays
  double bin_offset_; ///< m/z offset for the bins.
  bool use_flanks_; ///< use flanking peaks in calculation of XCorr theoretical

  /// used for xcorr
  double* observed_; ///< used for Xcorr: observed spectrum intensity array
  double* theoretical_; ///< used for Xcorr: theoretical spectrum intensity array

  /**
   * Initializes an empty scorer object
   */
  void init();

  /**
   * smooth all peaks in intensity array
   * Replaces the original array with the newly smooothed array
   */
  void smoothPeaks();

  /***
   * zero and extract peaks
   * extract peaks that are larger than mean + #step*stdev into new array
   * zero out the peaks that have been extracted
   * yes, the facter that a peak has removed will effect the fallowing peaks
   */
  void zeroPeakMeanStdev(
    double* original_array, ///< the array to normalize -in/out
    double* new_array, ///< the array to normalize -in/out                          
    int step                ///< is this 1 or 2 step -in
    );

  /**
   * \brief Zero and extract peaks
   *
   * Extract peaks that are larger than mean + #step*stdev into new
   * array.  Zero out the peaks that have been extracted.  Repeat twice,
   * than replace old array with extracted peak array.
   */
  void zeroPeaks();

  /**
   * keep only the peaks up to top rank peaks remove other peaks.
   * do second normalization on the top peaks back to max 100 intensity
   * replace old array with normalized top peak array
   */
  void extractPeaks(
    int top_rank  ///< keep the top ranking peaks -in
    );

  /**
   * equalize all peaks in a continous region to the largest peak within the continous bins
   * start from left to right
   */
  void equalizePeaks();

  /**
   * calculates all the necessay values for Sp score, related to the specfic ion_type
   * adds to intensity_sum and repeat_count
   *\returns the number of matches found from the predicted ions
   */
  int calculateIonTypeSp(
    IonSeries* ion_series, ///< the ion series to score against the spectrum -in
    double* intensity_sum,     ///< the total intensity sum of all matches so far -out
    ION_TYPE_T ion_type,      ///< the ion type to check -in
    int* repeat_count         ///< the repeated count of ions (ex. consecutive b ions) -out
    );

  /**
   * create the intensity array
   * SCORER must have been created for SP type
   * \returns true if successful, else FLASE
   */
  bool createIntensityArraySp(
    Crux::Spectrum* spectrum,    ///< the spectrum to score -in
    int charge               ///< the peptide charge -in 
    );

  /**
   * given a spectrum and ion series calculates the Sp score
   *\returns the sp score 
   */
  double genScoreSp(
    Crux::Spectrum* spectrum,    ///< the spectrum to score -in
    IonSeries* ion_series ///< the ion series to score against the spectrum -in
    );

  /*****************************************************
   * Xcorr related functions
   * 
   *****************************************************/

  /**
   * normalize each 10 regions of the observed spectrum to max 50
   */
  void normalizeEachRegion(
    double* observed,  ///< intensities to normalize
    double max_intensity_overall, /// the max intensity over entire spectrum
    double* max_intensity_per_region, ///< the max intensity in each 10 regions -in
    int region_selector ///< the size of each regions -in
    );

  /**
   * given a spectrum and ion series calculates the xcorr score
   *\returns the xcorr score 
   */
  double genScoreXcorr(
    Crux::Spectrum* spectrum,    ///< the spectrum to score -in
    IonSeries* ion_series ///< the ion series to score against the spectrum -in
    );

  /**
   * Create the intensity arrays for theoretical spectrum.
   * SCORER must have been created for XCORR type.
   * \returns true if successful, else FLASE
   */
  bool createIntensityArrayTheoretical(
    IonSeries* ion_series, ///< the ion series to score against the spectrum (theoretical) -in
    double*      theoretical ///< the empty theoretical spectrum -out
    );

  /*****************************************************
   * General purpose functions
   * 
   *****************************************************/
  /**
   * Creates the an array of ion constraints for GMTK models.
  
   * TODO do we need one for paired and single? Do we want an iterator?
   */
  IonConstraint** pairedIonConstraints();

  /**
   * Frees the paired ion_constraints array
   */
  void freePairedIonConstraints(
    IonConstraint** ion_constraints
    );

 public:

  /**
   * \returns An (empty) scorer object.
   */
  Scorer();

  /**
   * Instantiates a new scorer object from a SCORER_TYPE_T. 
   * \returns a new scorer object
   */
  Scorer(
    SCORER_TYPE_T type ///< the type of scorer -in
    );

  /**
   * Frees an allocated scorer object.
   */
  ~Scorer();

  /**
   * Score a spectrum vs. an ion series
   */
  double scoreSpectrumVIonSeries(
    Crux::Spectrum* spectrum,      ///< the spectrum to score -in
    IonSeries* ion_series ///< the ion series to score against the spectrum -in
  );

  /**
   * Frees the single_ion_constraints array
   */
  void freeSingleIonConstraints(
    IonConstraint** ion_constraints
    );

  /**
   * Creates the an array of ion constraints for GMTK models.
   * TODO do we need one for paired and single? Do we want an iterator?
   */
  IonConstraint** singleIonConstraints();

  /**
   * Score a spectrum vs. another spectrum
   */
  double scoreSpectrumVSpectrum(
    Crux::Spectrum* first_spectrum,   ///< the first spectrum to score -in
    Crux::Spectrum* second_spectrum   ///<  the second spectrum to score -in
  );

  /*******************************
   * get, set methods for scorer
   *******************************/

  /**
   *\returns the score type of the scorer
   */
  SCORER_TYPE_T getType();

  /**
   *sets the scorer type
   */
  void setType(
    SCORER_TYPE_T type ///< The type of scorer -in
    );

  /**
   *\returns the beta value of the scorer
   */
  double getSpBeta();

  /**
   *sets the scorer beta value
   */
  void setSpBeta(
    double sp_beta ///< used for Sp: the beta variable -in
    );

  /**
   *\returns the gamma value of the scorer
   */
  double getSpGamma();

  /**
   *sets the scorer gamma value
   */
  void setSpGamma(
    double sp_gamma ///< used for Sp: the gamma variable -in
    );


  /**
   *\returns the min_mz value of the scorer
   */
  double getSpMinMz();

  /**
   *sets the scorer min_mz value
   */
  void setSpMinMz(
    double sp_min_mz ///< used for Sp: the min_mz variable -in
    );


  /**
   *\returns the max_mz value of the scorer
   */
  double getSpMaxMz();

  /**
   *sets the scorer max_mz value
   */
  void setSpMaxMz(
    double sp_max_mz ///< used for Sp: the max_mz variable -in
    );

  /**
   *\returns the max bin index of the scorer array(s).
   */
  int getMaxBin();

  /**
   *\returns the sp_array_resolution value of the scorer
   */
  double getSpArrayResolution();

  /**
   *sets the scorer sp_array_resolution value
   */
  void setSpArrayResolution(
    double sp_array_resolution ///< used for Sp: the sp_array_resolution variable -in
    );

  /**
   *\returns the sp_sum_resolution value of the scorer
   */
  double getSpSumResolution();

  /**
   *sets the scorer sp_sum_resolution value
   */
  void setSpSumResolution(
    double sp_sum_resolution ///< used for Sp: the sp_sum_resolution variable -in
    );

  /**
   *\returns the equalize_resolution value of the scorer
   */
  double getSpEqualizeResolution();

  /**
   *sets the scorer equalize_resolution value
   */
  void setSpEqualizeResolution(
    double sp_equalize_resolution ///< used for Sp: the equalize_resolution variable -in
    );

  /**
   *\returns the fraction of b,y ions matched for scoring SP, the values is valid for the last ion series scored with this scorer object
   */
  double getSpBYIonFractionMatched();

  /**
   *\returns the number of possible matched b,y ions for scoring SP
   */
  int getSpBYIonMatched();

  /**
   *\returns the number of matched b,y ions for scoring SP
   */
  int getSpBYIonPossible();

  /**
   * Generate the processed peaks for the spectrum and return via the
   * intensities array.  It's implemented here so that
   * create_intensity_array_observed() can remain private and so that
   * the scorer->observed array can be accessed directly.
   * .
   */
  static void getProcessedPeaks(
    Crux::Spectrum* spectrum, 
    int charge,
    SCORER_TYPE_T score_type,
    double** intensities, ///< pointer to array of intensities
    int* mz_bins);

  /**
   * create the intensity arrays for both observed and theoretical spectrum
   * SCORER must have been created for XCORR type
   * \returns true if successful, else false
   */
  bool createIntensityArrayXcorr(
    Crux::Spectrum* spectrum,    ///< the spectrum to score(observed) -in
    int charge               ///< the peptide charge -in 
    );

  /**
   * Uses an iterative cross correlation
   *
   *\return the final cross correlation score between the observed and the
   *theoretical spectra
   */
  double crossCorrelation(
    double* theoretical ///< the theoretical spectrum to score against the observed spectrum -in
    );

  double* getIntensityArrayObserved();

  bool createIntensityArrayObserved(
    Crux::Spectrum* spectrum,    ///< the spectrum to score(observed) -in
    int charge               ///< the peptide charge -in 
    );

  /**
   * adds the intensity at add_idx
   * if, there already exist a peak at the index, only overwrite if
   * intensity is larger than the existing peak.
   */
  static void addIntensity(
    double* intensity_array, ///< the intensity array to add intensity at index add_idx -out
    int add_idx,            ///< the idex to add the intensity -in
    double intensity         ///< the intensity to add -in
    );
};


/*************************************
 * Score for LOGP_*
 ************************************/

/**
 * Compute a p-value for a given score w.r.t. an exponential with the given parameters.
 *\returns the -log(p_value) of the exponential distribution
 */
double score_logp_exp_sp(
  double sp_score, ///< The sp score for the scoring peptide -in
  double mean      ///< The overall mean of the sp scored peptides -in
  );

/**
 * Compute a p-value for a given score w.r.t. an exponential with the given parameters.
 *\returns the -log(p_value) of the exponential distribution with Bonferroni correction
 */
double score_logp_bonf_exp_sp(
  double sp_score, ///< The sp score for the scoring peptide -in
  double mean,      ///< The overall mean of the sp scored peptides -in
  int num_peptide  ///< The number of peptides scored for sp
  );

/**
 * Apply a Bonferroni correction to a given p-value.
 * \returns the corrected -log(p_value)
 */
double bonferroni_correction(
  double p_value, ///< The uncorrected p-value.
  int num_tests ///< The number of tests performed.
  );

/**
 * Compute a p-value for a given score w.r.t. a Weibull with given parameters.
 *\returns the p_value
 */
double compute_weibull_pvalue(
  double score, ///< The score for the scoring peptide -in
  double eta,   ///< The eta parameter of the Weibull -in
  double beta,  ///< The beta parameter of the Weibull -in
  double shift  ///< The shift parameter of the Weibull -in
  );

/**
 * Compute a p-value for a given score w.r.t. a Weibull with given parameters.
 *\returns the -log(p_value)
 */
double score_logp_bonf_weibull(
  double score, ///< The score for the scoring peptide -in
  double eta,  ///< The eta parameter of the Weibull
  double beta, ///< The beta parameter of the Weibull
  double shift, ///< The shift parameter of the Weibull
  int num_peptides ///< The number of peptides
  );


/**
 * Compute a p-value for a given score w.r.t. an EVD with the given parameters.
 *\returns the -log(p_value) of the EVD distribution 
 */
double score_logp_evd_xcorr(
  double xcorr_score, ///< The xcorr score for the scoring peptide -in
  double mu, ///<  EVD parameter Xcorr(characteristic value of extreme value distribution) -in
  double l_value ///< EVD parameter Xcorr(decay constant of extreme value distribution) -in
  );

/**
 * Compute a p-value for a given score w.r.t. an EVD with the given parameters.
 *\returns the -log(p_value) of the EVD distribution with Bonferroni correction
 */
double score_logp_bonf_evd_xcorr(
  double xcorr_score, ///< The xcorr score for the scoring peptide -in
  double mu, ///<  EVD parameter Xcorr(characteristic value of extreme value distribution) -in
  double l_value, ///< EVD parameter Xcorr(decay constant of extreme value distribution) -in
  int num_peptide  ///< The number of peptides scored for sp -in
  );


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
