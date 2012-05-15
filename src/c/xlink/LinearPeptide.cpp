
#include "LinearPeptide.h"
#include "XLink.h"
#include "ModifiedPeptidesIterator.h"
#include "IonSeries.h"

#include <iostream>

using namespace std;

LinearPeptide::LinearPeptide() {
  peptide_ = NULL;
  sequence_ = NULL;
}

LinearPeptide::LinearPeptide(char* sequence) {
  peptide_ = NULL;
  sequence_ = sequence;
}

LinearPeptide::LinearPeptide(Peptide* peptide) {
  peptide_ = peptide;
  sequence_ = NULL;
}


LinearPeptide::~LinearPeptide() {
}

void LinearPeptide::addCandidates(
  FLOAT_T min_mass,
  FLOAT_T max_mass,
  Index* index, 
  Database* database,
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods,
  XLinkMatchCollection& candidates) {

  int max_missed_cleavages = get_int_parameter("missed-cleavages");

  for (int mod_idx=0;mod_idx<num_peptide_mods; mod_idx++) {
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];

    //
    ModifiedPeptidesIterator* peptide_iterator =
      new ModifiedPeptidesIterator(min_mass, max_mass, peptide_mod, 
        false, index, database);


    while (peptide_iterator->hasNext()) {
      Peptide* peptide = peptide_iterator->next();
      XLinkMatch* new_candidate = new LinearPeptide(peptide);
      if (new_candidate->getNumMissedCleavages() <= max_missed_cleavages) {
        //cerr <<"Adding Linear Peptide:"<<new_candidate -> getSequenceString()<<" "<<new_candidate->getMass()<<endl;
        candidates.add(new_candidate);
        XLink::addAllocatedPeptide(peptide);
      } else {
        delete new_candidate;
        delete peptide;
      }
    }
    delete peptide_iterator;
  }
}

XLINKMATCH_TYPE_T LinearPeptide::getCandidateType() {
  return LINEAR_CANDIDATE;
}

string LinearPeptide::getSequenceString() {
  ostringstream oss;
  if (peptide_ == NULL) {
    oss << sequence_;

  } else {
    char* seq = peptide_->getModifiedSequenceWithMasses(MOD_MASSES_SEPARATE);
    oss << seq;
    free(seq);
  }

  oss << " ()";
  return oss.str();

}

FLOAT_T LinearPeptide::calcMass(MASS_TYPE_T mass_type) {
  if (peptide_ == NULL) {
    return Peptide::calcSequenceMass(sequence_, mass_type);
  } else {
  //  carp(CARP_INFO,"Calculating modified peptide mass");


   //FLOAT_T ans1 = calc_modified_peptide_mass(peptide_, mass_type);
    
  //  char* seq = get_peptide_sequence(peptide_);

  //  FLOAT_T ans2 = calc_sequence_mass(seq,mass_type);
  //  free(seq);

 //   FLOAT_T ans3 = get_peptide_peptide_mass(peptide_);

 //   carp(CARP_INFO,"%s %f %f %f",sequence_,ans1,ans2,ans3);
    
    return peptide_->calcModifiedMass(mass_type);
  }
}

XLinkMatch* LinearPeptide::shuffle() {

  Peptide* decoy_peptide = new Peptide(peptide_);

  decoy_peptide->transformToDecoy();

  XLink::addAllocatedPeptide(decoy_peptide);

  LinearPeptide* decoy = new LinearPeptide(decoy_peptide);

  return (XLinkMatch*)decoy;
}

void LinearPeptide::predictIons(IonSeries* ion_series, int charge) {

  char* seq = NULL;
  MODIFIED_AA_T* mod_seq = NULL;
  if (peptide_ == NULL) {
    seq = my_copy_string(sequence_);
    convert_to_mod_aa_seq(seq, &mod_seq); 
  } else {
    seq = peptide_->getSequence();
    mod_seq = peptide_->getModifiedAASequence();
  }
  ion_series->setCharge(charge);
  ion_series->update(seq, mod_seq);
  ion_series->predictIons();
/*
  IonFilteredIterator* ion_iterator = new IonFilteredIterator(ion_series, ion_constraint);
  while(ion_iterator->hasNext()) {
  
    Ion* ion = ion_iterator->next();
    if (isnan(ion->getMassZ())) {
      carp(CARP_FATAL, "NAN6");
    }
  }

*/

  free(seq);
  free(mod_seq);

}

string LinearPeptide::getIonSequence(Ion* ion) {

  string seq_str = string(sequence_);

  int cleavage_idx = ion->getCleavageIdx();
  if (ion->isForwardType() == B_ION) {
    return seq_str.substr(0,cleavage_idx);
  } else {
    return seq_str.substr(seq_str.length()-cleavage_idx,seq_str.length());
  }
}

Peptide* LinearPeptide::getPeptide(int peptide_idx) {
  if (peptide_idx == 0) {
    return peptide_;
  } else {
    return NULL;
  }
}

int LinearPeptide::getNumMissedCleavages() {
  set<int> skip;
  return peptide_->getMissedCleavageSites(skip);
}

bool LinearPeptide::isModified() {

  return peptide_->isModified();
}
