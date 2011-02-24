
#include "LinearPeptide.h"
#include "XLink.h"
#include "modified_peptides_iterator.h"
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

LinearPeptide::LinearPeptide(PEPTIDE_T* peptide) {
  peptide_ = peptide;
  sequence_ = NULL;
}


LinearPeptide::~LinearPeptide() {
}

void LinearPeptide::addCandidates(
  FLOAT_T min_mass,
  FLOAT_T max_mass,
  INDEX_T* index, 
  DATABASE_T* database,
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods,
  MatchCandidateVector& candidates) {

  int max_missed_cleavages = get_int_parameter("max-missed-cleavages");

  for (int mod_idx=0;mod_idx<num_peptide_mods; mod_idx++) {
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];

    //
    MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator =
      new_modified_peptides_iterator_from_mass_range(
						   min_mass,
						   max_mass,
						   peptide_mod,
						   FALSE,
						   index, database);


    while (modified_peptides_iterator_has_next(peptide_iterator)) {
      PEPTIDE_T* peptide = modified_peptides_iterator_next(peptide_iterator);
      MatchCandidate* new_candidate = new LinearPeptide(peptide);
      if (new_candidate->getNumMissedCleavages() <= max_missed_cleavages) {
        //cerr <<"Adding Linear Peptide:"<<new_candidate -> getSequenceString()<<" "<<new_candidate->getMass()<<endl;
        candidates.add(new_candidate);
        XLink::addAllocatedPeptide(peptide);
      } else {
        delete new_candidate;
        free_peptide(peptide);
      }
    }
    free_modified_peptides_iterator(peptide_iterator);
  }
}

MATCHCANDIDATE_TYPE_T LinearPeptide::getCandidateType() {
  return LINEAR_CANDIDATE;
}

string LinearPeptide::getSequenceString() {
  if (peptide_ == NULL) {
    string svalue(sequence_);
    return svalue;
  } else {
    char* seq = get_peptide_modified_sequence_with_masses(peptide_, FALSE);
    string svalue(seq);
    free(seq);
    return svalue;
  }
}

FLOAT_T LinearPeptide::calcMass(MASS_TYPE_T mass_type) {
  if (peptide_ == NULL) {
    return calc_sequence_mass(sequence_,mass_type);
  } else {
    return calc_modified_peptide_mass(peptide_, mass_type);
  }
}

MatchCandidate* LinearPeptide::shuffle() {

  PEPTIDE_T* decoy_peptide = copy_peptide(peptide_);

  transform_peptide_to_decoy(decoy_peptide);

  XLink::addAllocatedPeptide(decoy_peptide);

  LinearPeptide* decoy = new LinearPeptide(decoy_peptide);

  return (MatchCandidate*)decoy;
}

void LinearPeptide::predictIons(IonSeries* ion_series, int charge) {

  char* seq = NULL;
  MODIFIED_AA_T* mod_seq = NULL;
  if (peptide_ == NULL) {
    seq = sequence_;
    mod_seq = NULL; 
  } else {
    seq = get_peptide_sequence(peptide_);
    mod_seq = get_peptide_modified_aa_sequence(peptide_);
  }
  ion_series->setCharge(charge);
  ion_series->update(seq, mod_seq);
  ion_series->predictIons();
  free(seq);
  free(mod_seq);

}

string LinearPeptide::getIonSequence(Ion* ion) {
  return ion->getPeptideSequence();
}

PEPTIDE_T* LinearPeptide::getPeptide(int peptide_idx) {
  if (peptide_idx == 0) {
    return peptide_;
  } else {
    return NULL;
  }
}

int LinearPeptide::getNumMissedCleavages() {
  set<int> skip;
  return get_peptide_missed_cleavage_sites(peptide_, skip);
}

bool LinearPeptide::isModified() {

  return get_peptide_is_modified(peptide_);
}
