#include "MatchCandidateVector.h"
#include "XLinkPeptide.h"
#include "LinearPeptide.h"
#include "SelfLoopPeptide.h"
#include "XLinkScorer.h"

#include "Spectrum.h"

#include <iostream>


static const FLOAT_T MIN_XCORR_SHIFT = -5.0;
static const FLOAT_T MAX_XCORR_SHIFT  = 5.0;
//#define CORR_THRESHOLD 0.995   // Must achieve this correlation, else punt.
static const FLOAT_T CORR_THRESHOLD = 0.0;       // For now, turn off the threshold.
static const FLOAT_T XCORR_SHIFT = 0.05;

using namespace std;

void get_min_max_mass(
  FLOAT_T precursor_mz, 
  SpectrumZState& zstate, 
  FLOAT_T window,
  WINDOW_TYPE_T precursor_window_type,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass) {

  //cerr <<"mz: "
  //     <<precursor_mz
  //     <<" charge:"
  //     <<charge
  //     <<" mass:"<<mass
  //     <<" window:"<<window<<endl;
  if (precursor_window_type == WINDOW_MASS) {
    //cerr<<"WINDOW_MASS"<<endl;
    min_mass = zstate.getNeutralMass() - window;
    max_mass = zstate.getNeutralMass() + window;
  } else if (precursor_window_type == WINDOW_MZ) {
    //cerr<<"WINDOW_MZ"<<endl;
    double min_mz = precursor_mz - window;
    double max_mz = precursor_mz + window;
    min_mass = (min_mz - MASS_PROTON) * (double)zstate.getCharge();
    max_mass = (max_mz - MASS_PROTON) * (double)zstate.getCharge();
  } else if (precursor_window_type == WINDOW_PPM) {
    //cerr<<"WINDOW_PPM"<<endl;
    min_mass = zstate.getNeutralMass() / (1.0 + window * 1e-6);
    max_mass = zstate.getNeutralMass() / (1.0 - window * 1e-6);
  }
  
  //cerr<<"min:"<<min_mass<<" "<<"max: "<<max_mass<<endl;

}

void get_min_max_mass(
  FLOAT_T precursor_mz, 
  SpectrumZState& zstate,
  BOOLEAN_T use_decoy_window,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass) {
  
  if (use_decoy_window) {
    get_min_max_mass(precursor_mz,
		     zstate,
		     get_double_parameter("precursor-window-decoy"),
		     get_window_type_parameter("precursor-window-type-decoy"),
		     min_mass,
		     max_mass);
  } else {
    get_min_max_mass(precursor_mz,
		     zstate,
		     get_double_parameter("precursor-window"),
		     get_window_type_parameter("precursor-window-type"),
		     min_mass,
		     max_mass);
  }
}





bool compareXCorr(MatchCandidate* mc1, MatchCandidate* mc2) {
  return mc1->getXCorr() > mc2->getXCorr(); 
}

bool compareSP(MatchCandidate* mc1, MatchCandidate* mc2) {
  return mc1->getSP() > mc2->getSP();
}

MatchCandidateVector::MatchCandidateVector() {
  scan_ = 0;
}

MatchCandidateVector::MatchCandidateVector(MatchCandidateVector& vector) : std::vector<MatchCandidate*>() {
  

  precursor_mz_ = vector.precursor_mz_;
  zstate_ = vector.zstate_;
  scan_ = vector.scan_;

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
  FLOAT_T precursor_mz, 
  SpectrumZState& zstate,
  XLinkBondMap& bondmap,
  INDEX_T* index, 
  DATABASE_T* database, 
  PEPTIDE_MOD_T** peptide_mods, 
  int num_peptide_mods, 
  BOOLEAN_T use_decoy_window) {

  precursor_mz_ = precursor_mz;

  zstate_ = zstate;

  FLOAT_T min_mass;
  FLOAT_T max_mass;

  get_min_max_mass(precursor_mz, zstate, use_decoy_window, min_mass, max_mass);

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
  decoy_vector.precursor_mz_ = precursor_mz_;
  decoy_vector.zstate_ = zstate_;
  decoy_vector.scan_ = scan_;
  for (unsigned int idx=0;idx<size();idx++) {
    decoy_vector.add(at(idx)->shuffle());
  }
}

void MatchCandidateVector::scoreSpectrum(Spectrum* spectrum) {

  int max_ion_charge = get_max_ion_charge_parameter("max-ion-charge");

  carp(CARP_DEBUG, "Creating scorer");
  XLinkScorer scorer(
    spectrum, 
    min(zstate_.getCharge(), max_ion_charge));

  for (unsigned int idx=0;idx<size();idx++) {
    carp(CARP_DEBUG, "Scoring candidate:%d", idx);
    scorer.scoreCandidate(at(idx));
  }
  carp(CARP_DEBUG, "Done scoreSpectrum");
}

void MatchCandidateVector::sortByXCorr() {
  sort(begin(), end(), compareXCorr);
}

void MatchCandidateVector::sortBySP() {
  sort(begin(), end(), compareSP);
}


void MatchCandidateVector::setRanks() {
  sortBySP();

  int current_rank = 1;
  FLOAT_T last_score = at(0)->getSP();
  at(0)->setSPRank(current_rank);

  for (unsigned int candidate_idx = 1; candidate_idx < size() ; candidate_idx++) {
    FLOAT_T current_score = at(candidate_idx)->getSP();
    if (last_score != current_score) {
      current_rank++;
      last_score = current_score;
    }
    at(candidate_idx)->setSPRank(current_rank);
  }

  sortByXCorr();

  current_rank = 1;
  last_score = at(0)->getXCorr();
  at(0)->setXCorrRank(current_rank);
  for (unsigned int candidate_idx = 1; candidate_idx < size(); candidate_idx++) {
    FLOAT_T current_score = at(candidate_idx)->getXCorr();
    if (last_score != current_score) {
      current_rank++;
      last_score = current_score;
    }
    at(candidate_idx)->setXCorrRank(current_rank);
  }
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
  return zstate_.getCharge();
}

FLOAT_T MatchCandidateVector::getPrecursorMZ() {
  return precursor_mz_;
}

FLOAT_T MatchCandidateVector::getSpectrumNeutralMass() {
  return zstate_.getNeutralMass();
}


