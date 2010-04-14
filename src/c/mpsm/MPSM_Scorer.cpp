#include "MPSM_Scorer.h"
#include "scorer.h"
#include "modifications.h"

#include <iostream>
using namespace std;

MPSM_Scorer::MPSM_Scorer(SPECTRUM_T* spectrum, int max_charge) {
  spectrum_ = spectrum;
  max_charge_ = max_charge;

  scorer_ = new_scorer(XCORR);

  if (!create_intensity_array_xcorr(spectrum_, scorer_, max_charge_)) {
    carp(CARP_FATAL," Failed to produce observed intensity");
  }
  max_mz_ = (int)get_scorer_sp_max_mz(scorer_);
  cout <<"Max mz_:"<<max_mz_<<endl;
  theoretical_ = (FLOAT_T*)mymalloc(max_mz_* sizeof(FLOAT_T));

  ion_constraints_ = new ION_CONSTRAINT_T*[max_charge_];

  for (int charge =0; charge < max_charge_; charge++) {
     ion_constraints_[charge] = new_ion_constraint_smart(XCORR, charge+1);  
  }

  
}

MPSM_Scorer::~MPSM_Scorer() {
  free(theoretical_);
  free_scorer(scorer_);

  for (int charge=0;charge<max_charge_;charge++) {
    free_ion_constraint(ion_constraints_[charge]);
  }

  delete []ion_constraints_;

}

MPSM_Scorer* global_scorer = NULL;
SPECTRUM_T* global_spectrum = NULL;
int global_charge = 0;

FLOAT_T MPSM_Scorer::scoreMPSM(MPSM_Match& mpsm_match, SCORER_TYPE_T match_mode) {
  //only allow XCorr to be scored right now.

  if (match_mode != XCORR) {
    carp(CARP_FATAL,"Score type not implemented for mpsms!");
  }

  if (mpsm_match.numMatches() == 1) {
    //this should have already been calculated.
    MATCH_T* match = mpsm_match.getMatch(0);
    FLOAT_T score = get_match_score(match, match_mode);
    return score;
  } else {
    //calculate the joint.
    //assume the spectrum is the same for all matches.
    //can test with an assertion later.
    MATCH_T* match = mpsm_match.getMatch(0);
    SPECTRUM_T* spectrum = get_match_spectrum(match);
    int max_charge = mpsm_match.getChargeIndex().max();

    if (spectrum != global_spectrum || global_charge != max_charge) {
      cout <<"Creating new scorer:"<<endl;
      if (global_scorer != NULL) {
        delete global_scorer;
      }
      global_scorer = new MPSM_Scorer(spectrum, max_charge);
      global_spectrum = spectrum;
      global_charge = max_charge;
    }
    //cout <<"Scoring mpsm"<<endl;
    return global_scorer -> calcScore(mpsm_match, match_mode);
  }
}

FLOAT_T MPSM_Scorer::calcScore(MPSM_Match& mpsm_match, SCORER_TYPE_T match_mode) {
  //build the theoretical.
  //cout <<"Memset:"<<max_mz_<<endl;
  memset(theoretical_, 0, max_mz_ * sizeof(FLOAT_T));

  for (int match_idx = 0;
    match_idx < mpsm_match.numMatches();
    match_idx++) {
    //cout <<"Calculating theoretical for match:"<<match_idx<<endl;
    MATCH_T* current_match = mpsm_match.getMatch(match_idx);
    int current_charge = get_match_charge(current_match);
    char* sequence = get_match_sequence(current_match);
    MODIFIED_AA_T* mod_seq = get_match_mod_sequence(current_match);

    //cout <<"current charge:"<<current_charge<<endl;
    //cout <<"sequence:"<<sequence<<endl;

    ION_SERIES_T* ion_series = new_ion_series_generic(ion_constraints_[current_charge-1], 
      current_charge);
    update_ion_series(ion_series, sequence, mod_seq);
    predict_ions(ion_series);
    //count on this working, by just adding more ions to the theoretical
    //spectrum.
   //cout <<"Calling create_intensity_array_theoretical"<<endl;
   create_intensity_array_theoretical(scorer_, ion_series, theoretical_);
    //cout <<"Freeing objects"<<endl;
    free_ion_series(ion_series);
    free(sequence);
    free(mod_seq);
  }
  //cout <<"Calculating cross_correlation"<<endl;
  FLOAT_T score = cross_correlation(scorer_, theoretical_);
  mpsm_match.setScore(XCORR, score);     

  return score;
}


