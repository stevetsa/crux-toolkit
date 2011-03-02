#include "MPSM_Scorer.h"
#include "scorer.h"
#include "modifications.h"

#include <iostream>
using namespace std;



MPSM_Scorer::MPSM_Scorer() {
/*
  spectrum_ = NULL;
  max_charge_ = 0;
  scorer_ = NULL;
  ion_series_ = NULL;
  current_mpsm_ = NULL;
*/
}



FLOAT_T MPSM_Scorer::score(
  MPSM_Match& mpsm_match,
  SCORER_TYPE_T match_mode) {

  //TODO - do some caching to speed this process up.
  SCORER_T* scorer = new_scorer(match_mode);

  Spectrum* spectrum = mpsm_match.getSpectrum();
  
  int max_charge = mpsm_match.getMaxCharge();


  IonConstraint* ion_constraint = 
    IonConstraint::newIonConstraintSmart(match_mode, max_charge);
    
  IonSeries* ion_series = 
    new IonSeries(ion_constraint, max_charge);

  createIonSeries(mpsm_match, ion_series);

  FLOAT_T score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);


  //clean up.
  delete ion_series;
  delete ion_constraint;
  free_scorer(scorer);


  return score;
}

void MPSM_Scorer::createIonSeries(
  vector<string>& sequences, 
  vector<int>& charges,
  IonSeries* ion_series) {


  IonConstraint* ion_constraint = ion_series->getIonConstraint();

  for (int idx=0;idx<sequences.size();idx++) {

    int charge = charges.at(idx);
    const char* seq = sequences.at(idx).c_str();
    MODIFIED_AA_T* mod_seq = NULL;
    convert_to_mod_aa_seq((char*)seq, &mod_seq);

    IonSeries* temp_ions = new IonSeries(ion_constraint, charge);
    temp_ions->update((char*)seq, mod_seq);
    temp_ions->predictIons();

    for (IonIterator iter = temp_ions->begin();
      iter != temp_ions->end();
      ++iter) {
      ion_series->addIon(*iter);
    }
    //free the temp_ion_series, but don't delete them.
    IonSeries::freeIonSeries(temp_ions, false);
    free(mod_seq);
  }

  ion_series->setIsPredicted(true);

}


void MPSM_Scorer::createIonSeries(MPSM_Match& mpsm_match, IonSeries* ion_series) {

  IonConstraint* ion_constraint = ion_series->getIonConstraint();

  for (int idx=0;idx<mpsm_match.numMatches();idx++) {

    int charge = mpsm_match.getCharge(idx);
    char* seq = mpsm_match.getSequence(idx);
    MODIFIED_AA_T* mod_seq = mpsm_match.getModifiedSequence(idx);

    IonSeries* temp_ions = new IonSeries(ion_constraint, charge);
    temp_ions->update(seq, mod_seq);
    temp_ions->predictIons();

    for (IonIterator iter = temp_ions->begin();
      iter != temp_ions->end();
      ++iter) {
      ion_series->addIon(*iter);
    }
    //free the temp_ion_series, but don't delete them.
    IonSeries::freeIonSeries(temp_ions, false);

    free(seq);
    free(mod_seq);

  }

}

MPSM_Scorer::~MPSM_Scorer() {
}
