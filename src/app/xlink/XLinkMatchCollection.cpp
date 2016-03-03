/**
 * \file XLinkMatchCollection.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Collection of possible xlink products
 *****************************************************************************/

#include "XLinkMatchCollection.h"
#include "XLinkPeptide.h"
#include "LinearPeptide.h"
#include "SelfLoopPeptide.h"
#include "XLinkScorer.h"

#include "model/Spectrum.h"
#include "util/GlobalParams.h"

#include <iostream>


static const FLOAT_T MIN_XCORR_SHIFT = -5.0;
static const FLOAT_T MAX_XCORR_SHIFT  = 5.0;
//#define CORR_THRESHOLD 0.995   // Must achieve this correlation, else punt.
static const FLOAT_T CORR_THRESHOLD = 0.5;
static const FLOAT_T XCORR_SHIFT = 0.05;

using namespace std;

void get_min_max_mass(
  FLOAT_T precursor_mz, 
  SpectrumZState& zstate, 
  int isotope,
  FLOAT_T window,
  WINDOW_TYPE_T precursor_window_type,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass,
  FLOAT_T& precursor_mass) {

  //cerr <<"mz: "
  //     <<precursor_mz
  //     <<" charge:"
  //     <<charge
  //     <<" mass:"<<mass
  //     <<" window:"<<window<<endl;

  precursor_mass = zstate.getNeutralMass() + (double)isotope*MASS_NEUTRON;

  if (precursor_window_type == WINDOW_MASS) {
    //cerr<<"WINDOW_MASS"<<endl;
    min_mass = precursor_mass - window;
    max_mass = precursor_mass + window;
  } else if (precursor_window_type == WINDOW_MZ) {
    //cerr<<"WINDOW_MZ"<<endl;
    double min_mz = precursor_mz - window;
    double max_mz = precursor_mz + window;
    min_mass = (min_mz - MASS_PROTON) * (double)zstate.getCharge();
    max_mass = (max_mz - MASS_PROTON) * (double)zstate.getCharge();
  } else if (precursor_window_type == WINDOW_PPM) {
    //cerr<<"WINDOW_PPM"<<endl;
    min_mass = precursor_mass * (1.0 - window * 1e-6);
    max_mass = precursor_mass * (1.0 + window * 1e-6);
  }
  
  //cerr<<"min:"<<min_mass<<" "<<"max: "<<max_mass<<endl;

}

void get_min_max_mass(
  FLOAT_T precursor_mz, 
  SpectrumZState& zstate,
  int isotope,
  bool use_decoy_window,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass,
  FLOAT_T& precursor_mass) {
  
  if (use_decoy_window) {
    get_min_max_mass(precursor_mz,
         zstate,
                     isotope,
         get_double_parameter("precursor-window-weibull"),
         string_to_window_type(get_string_parameter("precursor-window-type-weibull")),
         min_mass,
		     max_mass, precursor_mass);
  } else {
    get_min_max_mass(precursor_mz,
      zstate,
      isotope,
      GlobalParams::getPrecursorWindow(),
      GlobalParams::getPrecursorWindowType(),		     
      min_mass,
      max_mass, 
      precursor_mass);
  }
}

/**
 * Default constructor
 */
XLinkMatchCollection::XLinkMatchCollection() : MatchCollection() {
  carp(CARP_DEBUG, "XLinkMatchCollection():start");
  scan_ = 0;
}

/**
 * Copy constructor
 */
XLinkMatchCollection::XLinkMatchCollection(
  XLinkMatchCollection& vector ///<collection to copy
  ) : MatchCollection() {
  
  carp(CARP_DEBUG, "XLinkMatchCollection(XLinkMatchCollection):start");
  precursor_mz_ = vector.precursor_mz_;
  zstate_ = vector.zstate_;
  scan_ = vector.scan_;

  for (int idx=0;idx<vector.getMatchTotal();idx++) {
    XLinkMatch* currentCandidate = (XLinkMatch*)vector[idx];
    XLinkMatch* copyCandidate = NULL;
    switch (currentCandidate -> getCandidateType()) {
    case LINEAR_CANDIDATE:
    case DEADLINK_CANDIDATE:
      copyCandidate = 
  new LinearPeptide(*(LinearPeptide*)currentCandidate);
      break;
    case SELFLOOP_CANDIDATE:
      copyCandidate =
  new SelfLoopPeptide(*(SelfLoopPeptide*)currentCandidate);
      break;
    case XLINK_INTER_CANDIDATE:
    case XLINK_INTRA_CANDIDATE:
    case XLINK_INTER_INTRA_CANDIDATE:
      copyCandidate =
  new XLinkPeptide(*(XLinkPeptide*)currentCandidate);
      break;
    case INVALID_CANDIDATE:
      carp(CARP_ERROR, "Invalid candidate type.");
      exit(1);
    }
    add(copyCandidate);
  }


}

/**
 * Constructor that finds all possible candidates
 */
/*
XLinkMatchCollection::XLinkMatchCollection() {
  
  carp(CARP_DEBUG, "XLinkMatchCollection(...)");

  FLOAT_T min_mass = get_double_parameter("min-mass");
  FLOAT_T max_mass = get_double_parameter("max-mass");

  addCandidates(
		NULL,
		0,
		1,
		min_mass, 
		max_mass,
		false);

}
*/

/**
 * Constructor that finds all candidates within a mass range
 */
void XLinkMatchCollection::addCandidates(
  Crux::Spectrum *spectrum,
  FLOAT_T precursor_mass,
  int precursor_charge,
  FLOAT_T min_mass, ///< min mass
  FLOAT_T max_mass, ///< max mass
  bool decoy
  ) {

  //carp(CARP_INFO, "XLinkMatchCollection.addCandidates() start");

  include_linear_peptides_ = get_boolean_parameter("xlink-include-linears");
  include_self_loops_ = get_boolean_parameter("xlink-include-selfloops");

  if (GlobalParams::getXLinkIncludeInter() ||
      GlobalParams::getXLinkIncludeIntra() ||
      GlobalParams::getXLinkIncludeInterIntra()) {
  
    carp(CARP_DEBUG, "Adding xlink candidates");
    carp(CARP_DEBUG, "precursor:%g", precursor_mass);
    carp(CARP_DEBUG, "min:%g", min_mass);
    carp(CARP_DEBUG, "max:%g", max_mass);
    XLinkPeptide::addCandidates(
      spectrum,
      precursor_mass,
      precursor_charge,
      min_mass, 
      max_mass,
      decoy,
      *this);
  }
  if (include_linear_peptides_) {

    LinearPeptide::addCandidates(
      min_mass,
      max_mass,
      decoy,
      *this);

  }

  if (include_self_loops_) {
  
    SelfLoopPeptide::addCandidates(
      min_mass,
      max_mass,
      decoy,
      *this);
  }
}


/**
 * Constructor for finding all candidates within a mass range
 */
XLinkMatchCollection::XLinkMatchCollection(
  Crux::Spectrum *spectrum, ///< spectrum
  SpectrumZState& zstate, ///< z-state
  bool decoy,
  bool use_decoy_window  
  ) {

  carp(CARP_DEBUG, "Inside XLinkMatchCollection....");

  precursor_mz_ = spectrum->getPrecursorMz();
  setZState(zstate);  


  FLOAT_T min_mass;
  FLOAT_T max_mass;
  const vector<int>& isotopes = GlobalParams::getIsotopeWindows(); 
  for (int idx = 0; idx < isotopes.size();idx++) {
    FLOAT_T precursor_mass;
    get_min_max_mass(precursor_mz_, zstate, isotopes[idx], use_decoy_window, min_mass, max_mass, precursor_mass);
    carp(CARP_DEBUG, "isotope %i precursor: %g min:%g max:%g", isotopes[idx], precursor_mass, min_mass, max_mass);
    addCandidates(spectrum, precursor_mass, zstate.getCharge(), 
		  min_mass, max_mass, decoy);
  }
}

/**
 * adds a candidate to the list
 */
void XLinkMatchCollection::add(
  XLinkMatch* candidate, ///< candidate to add
  bool copy
  ) {

  candidate->setZState(zstate_);
  candidate->setParent(this);
  addMatch(candidate);
  if (!copy) {
    candidate->decrementPointerCount();
  }
  experiment_size_++;

}

/**
 * \returns a candidate from the list by index
 */
XLinkMatch* XLinkMatchCollection::at(
  int idx ///< index of the candidate
  ) {
  if (idx < 0 || idx >= getMatchTotal()) {
    carp(CARP_FATAL, "XLinkMatchCollection:index %d out of bounds (0,%d)",
      idx, getMatchTotal());
  }
  return (XLinkMatch*)match_[idx];
}

/**
 * \returns a candidate from the list by index
 */
XLinkMatch* XLinkMatchCollection::operator [](
  int idx ///< index of the candidate
  ) {
  return (XLinkMatch*)match_[idx];
}


/**
 * shuffles the candidates and places the results in a decoy collection
 */
void XLinkMatchCollection::shuffle(
  XLinkMatchCollection& decoy_vector ///< collection to add decoys to
  ) {
  
  decoy_vector.precursor_mz_ = precursor_mz_;
  decoy_vector.zstate_ = zstate_;
  decoy_vector.scan_ = scan_;

  for (int idx=0;idx<getMatchTotal();idx++) {
    //cerr << "shuffling:" << at(idx)->getSequenceString()<<endl;
    decoy_vector.add(at(idx)->shuffle());
  }

}

/**
 * scores all candidates against the spectrum
 */
void XLinkMatchCollection::scoreSpectrum(
  Crux::Spectrum* spectrum ///< spectrum to score against
  ) {

  int max_ion_charge = get_max_ion_charge_parameter("max-ion-charge");

  carp(CARP_DEBUG, "Creating scorer");
  XLinkScorer scorer(
    spectrum, 
    min(zstate_.getCharge(), max_ion_charge));

  for (int idx=0;idx<getMatchTotal();idx++) {
    //carp(CARP_DEBUG, "Scoring candidate:%d", idx);
    //cerr << "XLinkMatchCollection::Scoreing" << at(idx)->getSequenceString()<<endl;
    scorer.scoreCandidate(at(idx));
  }

  // set the match_collection as having been scored
  scored_type_[XCORR] = true;
  if (get_boolean_parameter("compute-sp")) {
    scored_type_[SP] = true;
  }

  carp(CARP_DEBUG, "Done scoreSpectrum");
}


static int xcorrs_arr_length = 10;
static FLOAT_T *xcorrs  = (FLOAT_T*)mycalloc(10,
                                             sizeof(FLOAT_T));


/**
 * fits a weibull to the collection
 */
void XLinkMatchCollection::fitWeibull() {

  //create the array of x's and 
  shift_=0;
  eta_=0;
  beta_=0;
  correlation_=0;
  
  int nmatches = getMatchTotal();
  
  if (nmatches > xcorrs_arr_length) {
    carp(CARP_INFO, "getting larger array %d > %d", nmatches, xcorrs_arr_length);
    free(xcorrs);
    xcorrs = extractScores(XCORR);
    xcorrs_arr_length = nmatches;
  } else {
    extractScores(XCORR, xcorrs);
  }
  
  //FLOAT_T* xcorrs = extractScores(XCORR);
  // reverse sort the scores

  std::sort(xcorrs, xcorrs + nmatches, greater<FLOAT_T>());

  int num_tail_samples = (int)(nmatches * GlobalParams::getFractionToFit());

  FLOAT_T max_xcorr_shift = xcorrs[0];
  FLOAT_T min_xcorr_shift = min(MIN_XCORR_SHIFT, xcorrs[nmatches-1]);
  
  
  if (nmatches > 29) {
    max_xcorr_shift = xcorrs[29];
  
    carp(CARP_DEBUG, "min:%g max:%g", min_xcorr_shift, max_xcorr_shift);
    carp(CARP_DEBUG, "xcorr[%d]=%g", nmatches-1, xcorrs[nmatches-1]);
    carp(CARP_DEBUG, "xcorr[%d]=%g", 0, xcorrs[0]);
    fit_three_parameter_weibull(
      xcorrs,
      num_tail_samples,
      getMatchTotal(),
      min_xcorr_shift,
      max_xcorr_shift,
      XCORR_SHIFT,
      0,
      &eta_,
      &beta_,
      &shift_,
      &correlation_);
  }
  
  if (correlation_ < CORR_THRESHOLD) {
    carp(CARP_WARNING, "Weibull fit failed corr:%g n:%d nt:%d min:%g max:%g mins:%g maxs:%g",
      correlation_,
      nmatches,
      num_tail_samples,
      xcorrs[nmatches-1],
      xcorrs[0],
      min_xcorr_shift,
      max_xcorr_shift
      );
  } else {
    carp(CARP_DEBUG, "Fit success!");
  }
  
  
  //free(xcorrs);

}

/**
 * computes the p-value for the candidate
 */
void XLinkMatchCollection::computeWeibullPValue(
  int idx ///< candidate
  ) {
  if (correlation_ >= CORR_THRESHOLD) {
    at(idx)->computeWeibullPvalue(shift_, eta_, beta_);
  } else {
    //We could probably get a better estimate if we use the ecdf.
    carp(CARP_DEBUG, "Weibull fit failed corr:%g, assigning 1", correlation_);
    at(idx)->setPValue(1.0);
  }
}

/**
 * sets the scan for the collection
 */
void XLinkMatchCollection::setScan(
  unsigned int scan ///< scan number to set
  ) {
  scan_ = scan;
}

/**
 *\returns the scan number
 */
unsigned int XLinkMatchCollection::getScan() {
  return scan_;
}

/**
 *\returns the charge state of the collection
 */
int XLinkMatchCollection::getCharge() {
  return zstate_.getCharge();
}

/**
 *\returns the precursor m/z
 */
FLOAT_T XLinkMatchCollection::getPrecursorMZ() {
  return precursor_mz_;
}

/**
 * \returns the neutral mass of the collection
 */
FLOAT_T XLinkMatchCollection::getSpectrumNeutralMass() {
  return zstate_.getNeutralMass();
}

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */

