#include "XLinkScorer.h"
#include "scorer.h"
#include "ion.h"
#include "ion_series.h"

#include <iostream>

using namespace std;

XLinkScorer::XLinkScorer() {
  //  scorer = new_scorer(XCORR);
}

XLinkScorer::XLinkScorer(SPECTRUM_T* spectrum, int charge) {
  scorer_ = new_scorer(XCORR);
  spectrum_ = spectrum;
  charge_ = charge;
  ion_constraint_ = 
    new_ion_constraint_smart(XCORR, charge_);
  ion_series_ = 
    new_ion_series_generic(ion_constraint_, charge_);
}

XLinkScorer::~XLinkScorer() {
  //carp(CARP_INFO, "XLinkScorer::~XLinkScorer()");
  free(ion_series_);
  free(ion_constraint_);
  free_scorer(scorer_);
  
}

FLOAT_T XLinkScorer::scoreCandidate(MatchCandidate* candidate) {

  candidate->predictIons(ion_series_, charge_);
  FLOAT_T xcorr = score_spectrum_v_ion_series(scorer_, spectrum_, ion_series_);
  candidate->setXCorr(xcorr);
  //cerr<<candidate->getSequenceString()<<" xcorr "<<xcorr<<endl;
  return xcorr;

}
