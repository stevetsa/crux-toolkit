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

bool compareXCorr(MatchCandidate* mc1, MatchCandidate* mc2) {
  return mc1->getXCorr() > mc2->getXCorr(); 
}

MatchCandidateVector::MatchCandidateVector() {
  charge_ = 0;
  scan_ = 0;
  
  include_self_loops=FALSE;
  include_linear_peptides=FALSE;
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
  FLOAT_T precursor_mz, int charge, XLinkBondMap& bondmap,
  INDEX_T* index, DATABASE_T* database, PEPTIDE_MOD_T** peptide_mods, 
  int num_peptide_mods, BOOLEAN_T use_decoy_window) {
  
  charge_ = charge;
  precursor_mz_ = precursor_mz;

  include_linear_peptides=FALSE;
  include_self_loops=FALSE;

  //cerr<<"MatchCandidateVector(): Adding xlink candidates"<<endl;
    XLinkPeptide::addCandidates(precursor_mz, 
				charge, 
				bondmap,
				index,
				database,
				peptide_mods,
				num_peptide_mods,
				*this,
				use_decoy_window);
    
    if (include_linear_peptides) {
      //cerr <<"MatchCandidateVector(): Adding linear candidates"<<endl;
      LinearPeptide::addCandidates(precursor_mz,
				   charge,
				   index,
				   database,
				   peptide_mods,
				   num_peptide_mods,
				   *this,
				   use_decoy_window);
    }
    if (include_self_loops) {
      //cerr <<"MatchCandidateVector():  Adding self loop candidates"<<endl;
      SelfLoopPeptide::addCandidates(precursor_mz,
				     charge,
				     bondmap,
				     index,
				     database,
				     peptide_mods,
				     num_peptide_mods,
				     *this,
				     use_decoy_window);
    }

    //cerr<<"MatchCandidateVector(): Done adding candidates"<<endl;
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
  XLinkScorer scorer(spectrum, charge_);
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

  for (int i=0;i<3;i++) {
    carp(CARP_INFO,"%i:%f",i,xcorrs[i]);
  }

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
  /*
  eta_ = eta;
  beta_ = beta;
  shift_ = shift;
  */
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


