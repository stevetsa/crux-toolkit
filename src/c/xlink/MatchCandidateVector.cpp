#include "MatchCandidateVector.h"
#include "XLinkPeptide.h"
#include "LinearPeptide.h"
#include "SelfLoopPeptide.h"
#include "XLinkScorer.h"

#include "spectrum.h"

#include <iostream>


static const FLOAT_T MIN_XCORR_SHIFT = -5.0;
static const FLOAT_T MAX_XCORR_SHIFT  = 5.0;
//#define CORR_THRESHOLD 0.995   // Must achieve this correlation, else punt.
static const FLOAT_T CORR_THRESHOLD = 0.0;       // For now, turn off the threshold.
static const FLOAT_T XCORR_SHIFT = 0.05;

using namespace std;

void get_min_max_mass(
  FLOAT_T precursor_mz, 
  int charge, 
  FLOAT_T window,
  WINDOW_TYPE_T precursor_window_type,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass) {

  double mass = (precursor_mz - MASS_PROTON) * (double)charge;
  //cerr <<"mz: "
  //     <<precursor_mz
  //     <<" charge:"
  //     <<charge
  //     <<" mass:"<<mass
  //     <<" window:"<<window<<endl;
  if (precursor_window_type == WINDOW_MASS) {
    //cerr<<"WINDOW_MASS"<<endl;
    min_mass = mass - window;
    max_mass = mass + window;
  } else if (precursor_window_type == WINDOW_MZ) {
    //cerr<<"WINDOW_MZ"<<endl;
    double min_mz = precursor_mz - window;
    double max_mz = precursor_mz + window;
    min_mass = (min_mz - MASS_PROTON) * (double)charge;
    max_mass = (max_mz - MASS_PROTON) * (double)charge;
  } else if (precursor_window_type == WINDOW_PPM) {
    //cerr<<"WINDOW_PPM"<<endl;
    min_mass = mass / (1.0 + window * 1e-6);
    max_mass = mass / (1.0 - window * 1e-6);
  }
  
  //cerr<<"min:"<<min_mass<<" "<<"max: "<<max_mass<<endl;

}

void get_min_max_mass(
  FLOAT_T precursor_mz, 
  int charge, 
  BOOLEAN_T use_decoy_window,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass) {
  
  if (use_decoy_window) {
    get_min_max_mass(precursor_mz,
		     charge,
		     get_double_parameter("precursor-window-decoy"),
		     get_window_type_parameter("precursor-window-type-decoy"),
		     min_mass,
		     max_mass);
  } else {
    get_min_max_mass(precursor_mz,
		     charge,
		     get_double_parameter("precursor-window"),
		     get_window_type_parameter("precursor-window-type"),
		     min_mass,
		     max_mass);
  }
}





bool compareXCorr(MatchCandidate* mc1, MatchCandidate* mc2) {
  return mc1->getXCorr() > mc2->getXCorr(); 
}

MatchCandidateVector::MatchCandidateVector() {
  charge_ = 0;
  scan_ = 0;
}

MatchCandidateVector::MatchCandidateVector(MatchCandidateVector& vector) : std::vector<MatchCandidate*>() {
  
  charge_ = vector.charge_;
  scan_ = vector.scan_;
  precursor_mz_ = vector.precursor_mz_;
  for (unsigned int idx=0;idx<vector.size();idx++) {
    MatchCandidate* currentCandidate = vector[idx];
    MatchCandidate* copyCandidate = NULL;
    switch (currentCandidate -> getCandidateType()) {
    case LINEAR_CANDIDATE:
      copyCandidate = 
	new LinearPeptide(*(LinearPeptide*)currentCandidate);
      break;
    case SELFLOOP_CANDIDATE:
      copyCandidate =
	new SelfLoopPeptide(*(SelfLoopPeptide*)currentCandidate);
      break;
    case XLINK_CANDIDATE:
      copyCandidate =
	new XLinkPeptide(*(XLinkPeptide*)currentCandidate);
      break;
    }
    add(copyCandidate);
  }

}

MatchCandidateVector::MatchCandidateVector(
  XLinkBondMap& bondmap,
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods,
  INDEX_T* index,
  DATABASE_T* database) {


  FLOAT_T min_mass = get_double_parameter("min-mass");
  FLOAT_T max_mass = get_double_parameter("max-mass");

  addCandidates(
    min_mass, 
    max_mass, 
    bondmap, 
    index, 
    database, 
    peptide_mods, 
    num_peptide_mods);
  
}


void MatchCandidateVector::addCandidates(
  FLOAT_T min_mass,
  FLOAT_T max_mass,
  XLinkBondMap& bondmap,
  INDEX_T* index,
  DATABASE_T* database,
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods) {


  include_linear_peptides = get_boolean_parameter("xlink-include-linears");
  include_self_loops = get_boolean_parameter("xlink-include-selfloops");

  XLinkPeptide::addCandidates(
    min_mass, 
    max_mass,
    bondmap,
    index,
    database,
    peptide_mods,
    num_peptide_mods,
    *this);

  if (include_linear_peptides) {

    LinearPeptide::addCandidates(
      min_mass,
      max_mass,
      index,
      database,
      peptide_mods,
      num_peptide_mods,
      *this);

  }

  if (include_self_loops) {

    SelfLoopPeptide::addCandidates(
      min_mass,
      max_mass,
      bondmap,
      index,
      database,
      peptide_mods,
      num_peptide_mods,
      *this);
  }
}


MatchCandidateVector::MatchCandidateVector(
  FLOAT_T precursor_mz, int charge, XLinkBondMap& bondmap,
  INDEX_T* index, DATABASE_T* database, PEPTIDE_MOD_T** peptide_mods, 
  int num_peptide_mods, BOOLEAN_T use_decoy_window) {
  
  charge_ = charge;
  precursor_mz_ = precursor_mz;


  FLOAT_T min_mass;
  FLOAT_T max_mass;

  get_min_max_mass(precursor_mz, charge, use_decoy_window, min_mass, max_mass);

  addCandidates(min_mass, max_mass, bondmap, index, database, peptide_mods, num_peptide_mods);
}


MatchCandidateVector::~MatchCandidateVector() {
  for (unsigned int idx=0;idx<size();idx++) {
    delete at(idx);
  }
  clear();
}

void MatchCandidateVector::add(MatchCandidate* candidate) {
  push_back(candidate);
  candidate->setParent(this);
}

void MatchCandidateVector::shuffle() {
}

void MatchCandidateVector::shuffle(MatchCandidateVector& decoy_vector) {
  decoy_vector.charge_ = charge_;
  decoy_vector.scan_ = scan_;
  for (unsigned int idx=0;idx<size();idx++) {
    decoy_vector.add(at(idx)->shuffle());
  }
}

void MatchCandidateVector::scoreSpectrum(SPECTRUM_T* spectrum) {

  int max_ion_charge = get_max_ion_charge_parameter("max-ion-charge");
  
  XLinkScorer scorer(spectrum, min(charge_, max_ion_charge));
  for (unsigned int idx=0;idx<size();idx++) {
    scorer.scoreCandidate(at(idx));
  }
}

void MatchCandidateVector::sortByXCorr() {
  sort(begin(), end(), compareXCorr);
}

void MatchCandidateVector::fitWeibull(
  FLOAT_T& shift, 
  FLOAT_T& eta, 
  FLOAT_T& beta, 
  FLOAT_T& corr) {

  //create the array of x's and 
  shift=0;
  eta=0;
  beta=0;
  corr=0;

  FLOAT_T* xcorrs = new FLOAT_T[size()];

  for (unsigned int idx=0;idx<size();idx++) {
    xcorrs[idx] = at(idx)->getXCorr();
  }

  sort(xcorrs, xcorrs+size(), greater<FLOAT_T>());

  double fraction_to_fit = get_double_parameter("fraction-top-scores-to-fit");
  int num_tail_samples = (int)(size() * fraction_to_fit);

  fit_three_parameter_weibull(xcorrs,
			      num_tail_samples,
			      size(),
			      MIN_XCORR_SHIFT,
			      MAX_XCORR_SHIFT,
			      XCORR_SHIFT,
			      CORR_THRESHOLD,
			      &eta,
			      &beta,
			      &shift,
			      &corr);

  free(xcorrs);
}


void MatchCandidateVector::setScan(unsigned int scan) {
  scan_ = scan;
}

unsigned int MatchCandidateVector::getScan() {
  return scan_;
}

int MatchCandidateVector::getCharge() {
  return charge_;
}

FLOAT_T MatchCandidateVector::getPrecursorMZ() {
  return precursor_mz_;
}

FLOAT_T MatchCandidateVector::getSpectrumNeutralMass() {
  return (precursor_mz_-MASS_PROTON)*(double)charge_;
}


