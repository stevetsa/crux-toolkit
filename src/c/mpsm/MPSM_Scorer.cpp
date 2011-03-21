#include "MPSM_Scorer.h"
#include "scorer.h"
#include "modifications.h"

#include <iostream>
using namespace std;

static vector<vector< IonConstraint*> > ion_constraints_;
static map<IonConstraint*, vector<IonSeries*> > ion_series_map_;


MPSM_Scorer::MPSM_Scorer() {
/*
  spectrum_ = NULL;
  max_charge_ = 0;
  scorer_ = NULL;
  ion_series_ = NULL;
  current_mpsm_ = NULL;
*/
}

MPSM_Scorer::MPSM_Scorer(
  Spectrum* spectrum,
  int max_charge,
  SCORER_TYPE_T scorer_type) {


  spectrum_ = spectrum;
  max_charge_ = max_charge;
  scorer_type_ = scorer_type;

  scorer_ = new_scorer(scorer_type_);
  ion_constraint_ = getIonConstraint(scorer_type_,max_charge_);

}

IonConstraint* MPSM_Scorer::getIonConstraint(
  SCORER_TYPE_T scorer_type,
  int charge) {

  while (ion_constraints_.size() <= (int)scorer_type) {

    vector<IonConstraint*> ion_constraint_vec_;
    ion_constraints_.push_back(ion_constraint_vec_);
    
  }

  vector<IonConstraint*>& scorer_ion_constraints = ion_constraints_[(int)scorer_type];

  while (scorer_ion_constraints.size() < charge) {
    IonConstraint* ion_constraint =
      IonConstraint::newIonConstraintSmart(scorer_type, scorer_ion_constraints.size()+1);
    
    scorer_ion_constraints.push_back(ion_constraint);
  }

  return scorer_ion_constraints.at(charge-1);
  

}

IonSeries* MPSM_Scorer::getIonSeries(
  IonConstraint* ion_constraint,
  int charge) {

  map<IonConstraint*, vector<IonSeries*> >::iterator find_iter;
  find_iter = ion_series_map_.find(ion_constraint);

  if (find_iter == ion_series_map_.end()) {
    vector<IonSeries*> ion_series_vec;
    ion_series_map_[ion_constraint] = ion_series_vec;
  }

  vector<IonSeries*>& ion_series_vec = ion_series_map_[ion_constraint];

  while (ion_series_vec.size() < charge) {
    IonSeries* ion_series = new
      IonSeries(ion_constraint, ion_series_vec.size()+1);
    ion_series_vec.push_back(ion_series);
  }
    
  return ion_series_vec.at(charge-1);

}



FLOAT_T MPSM_Scorer::calcScore(
  MPSM_Match& mpsm_match
  ) {

  IonSeries* ion_series =
    new IonSeries(ion_constraint_, max_charge_);

  createIonSeries(mpsm_match, ion_series);

  FLOAT_T score = score_spectrum_v_ion_series(scorer_, spectrum_, ion_series);

  delete ion_series;
  
  return score;

}


FLOAT_T MPSM_Scorer::score(
  MPSM_Match& mpsm_match,
  SCORER_TYPE_T match_mode) {

  //TODO - do some caching to speed this process up.
  SCORER_T* scorer = new_scorer(match_mode);

  Spectrum* spectrum = mpsm_match.getSpectrum();
  
  int max_charge = mpsm_match.getMaxCharge();


  

  IonConstraint* ion_constraint = getIonConstraint(match_mode, max_charge);
    
  IonSeries* ion_series = 
    new IonSeries(ion_constraint, max_charge);

  createIonSeries(mpsm_match, ion_series);

  FLOAT_T score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);


  //clean up.
  delete ion_series;
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

    IonSeries* temp_ions = getIonSeries(ion_constraint, charge);
    temp_ions->update((char*)seq, mod_seq);
    temp_ions->predictIons();

    for (IonIterator iter = temp_ions->begin();
      iter != temp_ions->end();
      ++iter) {
      ion_series->addIon(*iter);
    }
    //free the temp_ion_series, but don't delete them.
    temp_ions->clearIons(false);
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

    IonSeries* temp_ions = getIonSeries(ion_constraint, charge);
    temp_ions->update(seq, mod_seq);
    temp_ions->predictIons();

    for (IonIterator iter = temp_ions->begin();
      iter != temp_ions->end();
      ++iter) {
      ion_series->addIon(*iter);
    }
    //free the temp_ion_series, but don't delete them.
    temp_ions->clearIons(false);

    free(seq);
    free(mod_seq);

  }

}

MPSM_Scorer::~MPSM_Scorer() {
  spectrum_ = NULL;
  max_charge_ = 0;
  free_scorer(scorer_);
  ion_constraint_ = NULL;
}
