#include "XLinkMatchCollection.h"
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
  bool use_decoy_window,
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




/*
bool compareXCorr(XLinkMatch* mc1, XLinkMatch* mc2) {
  return mc1->getXCorr() > mc2->getXCorr(); 
}

bool compareSP(XLinkMatch* mc1, XLinkMatch* mc2) {
  return mc1->getSP() > mc2->getSP();
}
*/

XLinkMatchCollection::XLinkMatchCollection() : MatchCollection () {
  scan_ = 0;
}

XLinkMatchCollection::XLinkMatchCollection(
  XLinkMatchCollection& vector
  ) : MatchCollection() {
  

  (void)vector;


  precursor_mz_ = vector.precursor_mz_;
  zstate_ = vector.zstate_;
  scan_ = vector.scan_;

  for (int idx=0;idx<vector.getMatchTotal();idx++) {
    XLinkMatch* currentCandidate = (XLinkMatch*)vector.match_[idx];
    XLinkMatch* copyCandidate = NULL;
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

XLinkMatchCollection::XLinkMatchCollection(
  XLinkBondMap& bondmap,
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods,
  Index* index,
  Database* database) {


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


void XLinkMatchCollection::addCandidates(
  FLOAT_T min_mass,
  FLOAT_T max_mass,
  XLinkBondMap& bondmap,
  Index* index,
  Database* database,
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


XLinkMatchCollection::XLinkMatchCollection(
  FLOAT_T precursor_mz, 
  SpectrumZState& zstate,
  XLinkBondMap& bondmap,
  Index* index, 
  Database* database, 
  PEPTIDE_MOD_T** peptide_mods, 
  int num_peptide_mods, 
  bool use_decoy_window) {

  precursor_mz_ = precursor_mz;

  zstate_ = zstate;

  FLOAT_T min_mass;
  FLOAT_T max_mass;

  get_min_max_mass(precursor_mz, zstate, use_decoy_window, min_mass, max_mass);

  addCandidates(min_mass, max_mass, bondmap, index, database, peptide_mods, num_peptide_mods);
}


XLinkMatchCollection::~XLinkMatchCollection()  {
//TODO - make sure we dont have to do any special handling here.
}

void XLinkMatchCollection::add(XLinkMatch* candidate) {

  addMatch(candidate);
  candidate->setParent(this);

}

void XLinkMatchCollection::shuffle() {
}

void XLinkMatchCollection::shuffle(XLinkMatchCollection& decoy_vector) {
  decoy_vector.precursor_mz_ = precursor_mz_;
  decoy_vector.zstate_ = zstate_;
  decoy_vector.scan_ = scan_;

  for (int idx=0;idx<getMatchTotal();idx++) {
    XLinkMatch* current = (XLinkMatch*)match_[idx];
    decoy_vector.addMatch((Match*)current->shuffle());
  }

}

void XLinkMatchCollection::scoreSpectrum(Spectrum* spectrum) {

  int max_ion_charge = get_max_ion_charge_parameter("max-ion-charge");

  carp(CARP_DEBUG, "Creating scorer");
  XLinkScorer scorer(
    spectrum, 
    min(zstate_.getCharge(), max_ion_charge));

  for (int idx=0;idx<getMatchTotal();idx++) {
    carp(CARP_DEBUG, "Scoring candidate:%d", idx);
    scorer.scoreCandidate((XLinkMatch*)match_[idx]);
  }

  carp(CARP_DEBUG, "Done scoreSpectrum");
}


void XLinkMatchCollection::setRanks() {
/*
  int current_rank = 1;
  FLOAT_T last_score = 0;

  if (get_boolean_parameter("compute-sp")) {
    sortBySP();

    current_rank = 1;
    last_score = at(0)->getSP();
    at(0)->setSPRank(current_rank);

    for (unsigned int candidate_idx = 1; candidate_idx < size() ; candidate_idx++) {
      FLOAT_T current_score = 0;//at(candidate_idx)->getSP();
      if (last_score != current_score) {
        current_rank++;
        last_score = current_score;
      }
      at(candidate_idx)->setSPRank(current_rank);
    }
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
*/
}


void XLinkMatchCollection::fitWeibull() {

  //create the array of x's and 
  shift_=0;
  eta_=0;
  beta_=0;
  correlation_=0;

  FLOAT_T* xcorrs = extractScores(XCORR);

  double fraction_to_fit = get_double_parameter("fraction-top-scores-to-fit");
  int num_tail_samples = (int)(0/*size()*/ * fraction_to_fit);

  fit_three_parameter_weibull(xcorrs,
			      num_tail_samples,
			      getMatchTotal(),
			      MIN_XCORR_SHIFT,
			      MAX_XCORR_SHIFT,
			      XCORR_SHIFT,
			      CORR_THRESHOLD,
			      &eta_,
			      &beta_,
			      &shift_,
			      &correlation_);

  myfree(xcorrs);
}

void XLinkMatchCollection::computeWeibullPValue(
  int idx
  ) {

  ((XLinkMatch*)match_[idx])->computeWeibullPvalue(shift_, eta_, beta_);

}


void XLinkMatchCollection::setScan(unsigned int scan) {
  scan_ = scan;
}

unsigned int XLinkMatchCollection::getScan() {
  return scan_;
}

int XLinkMatchCollection::getCharge() {
  return zstate_.getCharge();
}

FLOAT_T XLinkMatchCollection::getPrecursorMZ() {
  return precursor_mz_;
}

FLOAT_T XLinkMatchCollection::getSpectrumNeutralMass() {
  return zstate_.getNeutralMass();
}


