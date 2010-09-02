
#include "LinearPeptide.h"
#include "XLink.h"
#include "modified_peptides_iterator.h"
#include "ion_series.h"

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

void LinearPeptide::addCandidates(FLOAT_T precursor_mz, int charge, 
			    INDEX_T* index, DATABASE_T* database,
			    PEPTIDE_MOD_T** peptide_mods,
			    int num_peptide_mods,
			    MatchCandidateVector& candidates,
			     BOOLEAN_T use_decoy_window) {

  int max_missed_cleavages = get_int_parameter("max-missed-cleavages");

  FLOAT_T min_mass;
  FLOAT_T max_mass;

  get_min_max_mass(precursor_mz, charge, use_decoy_window, min_mass, max_mass);

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

      if (get_peptide_missed_cleavage_sites(peptide) <= max_missed_cleavages) {

        MatchCandidate* new_candidate = new LinearPeptide(peptide);
        //cerr <<"Adding Linear Peptide:"<<new_candidate -> getSequenceString()<<" "<<new_candidate->getMass()<<endl;
        candidates.add(new_candidate);
        XLink::addAllocatedPeptide(peptide);
      } else {
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

FLOAT_T LinearPeptide::getMass() {
  if (peptide_ == NULL) {
    return calc_sequence_mass(sequence_,AVERAGE);
  } else {
    return get_peptide_peptide_mass(peptide_);
  }
}

MatchCandidate* LinearPeptide::shuffle() {

  PEPTIDE_T* decoy_peptide = copy_peptide(peptide_);

  transform_peptide_to_decoy(decoy_peptide);

  XLink::addAllocatedPeptide(decoy_peptide);

  LinearPeptide* decoy = new LinearPeptide(decoy_peptide);

  return (MatchCandidate*)decoy;
}

void LinearPeptide::predictIons(ION_SERIES_T* ion_series, int charge) {

  char* seq = NULL;
  MODIFIED_AA_T* mod_seq = NULL;
  if (peptide_ == NULL) {
    seq = sequence_;
    mod_seq = NULL; 
  } else {
    seq = get_peptide_sequence(peptide_);
    mod_seq = get_peptide_modified_aa_sequence(peptide_);
  }
  set_ion_series_charge(ion_series, charge);
  update_ion_series(ion_series, seq, mod_seq);
  predict_ions(ion_series);
  free(seq);
  free(mod_seq);

}

string LinearPeptide::getIonSequence(ION_T* ion) {
  return get_ion_peptide_sequence(ion);
}
