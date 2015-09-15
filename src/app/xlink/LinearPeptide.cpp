
#include "LinearPeptide.h"
#include "XLink.h"
#include "model/ModifiedPeptidesIterator.h"
#include "model/IonSeries.h"
#include "util/GlobalParams.h"


#include <iostream>

using namespace std;


vector<LinearPeptide> LinearPeptide::target_linear_peptides_;
vector<LinearPeptide> LinearPeptide::decoy_linear_peptides_;


/**
 * Default constructor
 */
LinearPeptide::LinearPeptide() {
  peptide_ = NULL;
  sequence_ = NULL;
}

/**
 * Constructor with sequence
 */
LinearPeptide::LinearPeptide(
  char* sequence ///< sequence string
  ) {
  peptide_ = NULL;
  sequence_ = sequence;
}

/**
 * Constructor from a Crux Peptide
 */
LinearPeptide::LinearPeptide(
  Crux::Peptide* peptide ///< peptide object
  ) {
  peptide_ = peptide;
  sequence_ = NULL;
}


vector<LinearPeptide>::iterator LinearPeptide::getLinearPeptidesBegin(
  Database* database,
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods,
  bool is_decoy
  ) {

  if (is_decoy) {
    if (decoy_linear_peptides_.empty()) {
      generateAllLinearPeptides(decoy_linear_peptides_, database, peptide_mods, num_peptide_mods, is_decoy);
    }
    return(decoy_linear_peptides_.begin());
  } else {
    if (target_linear_peptides_.empty()) {
      generateAllLinearPeptides(target_linear_peptides_, database, peptide_mods, num_peptide_mods, is_decoy);
    }
    return(target_linear_peptides_.begin());
  }
}

vector<LinearPeptide>::iterator LinearPeptide::getLinearPeptidesEnd(
  Database* database, 
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods,
  bool is_decoy
  ){

  if (is_decoy) {
    if (decoy_linear_peptides_.empty()) {
      generateAllLinearPeptides(decoy_linear_peptides_, database, peptide_mods, num_peptide_mods, is_decoy);
    }
    return(decoy_linear_peptides_.end());
  } else {
    if (target_linear_peptides_.empty()) {
      generateAllLinearPeptides(target_linear_peptides_, database, peptide_mods, num_peptide_mods, is_decoy);
    }
    return(target_linear_peptides_.end());
  }

}

void LinearPeptide::generateAllLinearPeptides(
  vector<LinearPeptide>& linear_peptides,
  Database* database,
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods,
  bool decoy
  ) {

  linear_peptides.empty();
  int max_missed_cleavages = GlobalParams::getMissedCleavages();
  for (int mod_idx=0;mod_idx<num_peptide_mods; mod_idx++) {
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];
    double delta_mass = peptide_mod_get_mass_change(peptide_mod);
    //
    ModifiedPeptidesIterator* peptide_iterator =
      new ModifiedPeptidesIterator(
        GlobalParams::getMinMass(), 
        GlobalParams::getMaxMass(), 
        peptide_mod, 
        false, 
        database,
        0);

    //add the targets
    while (peptide_iterator->hasNext()) {
      
      Crux::Peptide* peptide = peptide_iterator->next();
      LinearPeptide linear_peptide(peptide);
      linear_peptide.getMass(GlobalParams::getIsotopicMass());
      if (linear_peptide.getNumMissedCleavages() <= max_missed_cleavages) {
        linear_peptides.push_back(linear_peptide);
	//XLink::addAllocatedPeptide(peptide);
	} else {
        delete peptide;
       }
      }
    delete peptide_iterator;
  }
  //carp(CARP_INFO, "sorting");
  sort(linear_peptides.begin(), linear_peptides.end(), compareLinearPeptideMass);

  
  //carp(CARP_INFO, "printing");
  //for (size_t idx=0;idx < linear_peptides.size();idx++) {
  //  string sequence = linear_peptides[idx].getSequenceString();
  //  carp(CARP_INFO, "%d %f %s", idx, linear_peptides[idx].getMassConst(MONO), sequence.c_str());
  //}
  //carp(CARP_FATAL, "stopping here");
  //  carp(CARP_FATAL, "Stopping here");

}

/**
 *Add candidates to the XLinkMatchCollection that are linear
 */
void LinearPeptide::addCandidates(
  FLOAT_T min_mass,  ///< min mass
  FLOAT_T max_mass,  ///< max mass
  Database* database, ///< protein database
  PEPTIDE_MOD_T** peptide_mods, ///< modifications peptide can take
  int num_peptide_mods, ///< Number of possible peptide mods
  bool is_decoy, ///< generate decoy canidates
  XLinkMatchCollection& candidates ///< Vector of candidate -inout
  ) {

  if (is_decoy) {
    if (decoy_linear_peptides_.empty()) {
      generateAllLinearPeptides(decoy_linear_peptides_, database, peptide_mods, num_peptide_mods, is_decoy);
    }
    vector<LinearPeptide>::iterator iter = 
      lower_bound(decoy_linear_peptides_.begin(), decoy_linear_peptides_.end(), min_mass, compareLinearPeptideMassToFLOAT);
    while (iter != decoy_linear_peptides_.end() && iter -> getMass(GlobalParams::getIsotopicMass()) <= max_mass) {
      iter->incrementPointerCount();
      candidates.add(&(*iter));
    }
  } else {
    if (target_linear_peptides_.empty()) {
      generateAllLinearPeptides(target_linear_peptides_, database, peptide_mods, num_peptide_mods, is_decoy);
    }
    vector<LinearPeptide>::iterator iter = 
      lower_bound(target_linear_peptides_.begin(), target_linear_peptides_.end(), min_mass, compareLinearPeptideMassToFLOAT); 
    while(iter != target_linear_peptides_.end() && iter -> getMass(GlobalParams::getIsotopicMass()) <= max_mass) {
     
      iter->incrementPointerCount(); 
      candidates.add(&(*iter));
      ++iter;
    }
  }
}

/**
 * returns the candidate type, either a deadlink or a linear candidate
 */
XLINKMATCH_TYPE_T LinearPeptide::getCandidateType() {
  if (isModified()) {
    return DEADLINK_CANDIDATE;
  } else {
    return LINEAR_CANDIDATE;
  }
}

/**
 * \returns the sequence of the peptide in string format
 */
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

/**
 * \returns the mass of the peptide
 */
FLOAT_T LinearPeptide::calcMass(
  MASS_TYPE_T mass_type ///< MONO or AVERAGE
  ) {
  
  if (peptide_ == NULL) {
    return Crux::Peptide::calcSequenceMass(sequence_, mass_type);
  } else {
    return peptide_->calcModifiedMass(mass_type);
  }
}

/**
 *\returns a shuffled version of the peptide
 */
XLinkMatch* LinearPeptide::shuffle() {
  string seq = getSequenceString();
  carp(CARP_DEBUG, "LinearPeptide::shuffle %s", seq.c_str());
  Crux::Peptide* decoy_peptide = new Crux::Peptide(peptide_);

  decoy_peptide->transformToDecoy();

  XLink::addAllocatedPeptide(decoy_peptide);

  LinearPeptide* decoy = new LinearPeptide(decoy_peptide);

  return (XLinkMatch*)decoy;
}

/**
 * predicts the ions for this peptide
 */
void LinearPeptide::predictIons(
  IonSeries* ion_series, ///< ion series to place the ions
  int charge ///< charge state of the peptide
  ) {

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
  free(seq);
  free(mod_seq);

}

/**
 *\returns the ion sequence as a string
 */
string LinearPeptide::getIonSequence(
  Ion* ion ///< ion object
  ) {

  string seq_str = string(sequence_);

  int cleavage_idx = ion->getCleavageIdx();
  if (ion->isForwardType() == B_ION) {
    return seq_str.substr(0,cleavage_idx);
  } else {
    return seq_str.substr(seq_str.length()-cleavage_idx,seq_str.length());
  }
}

/**
 * \returns the peptide for this match
 */
Crux::Peptide* LinearPeptide::getPeptide(int peptide_idx) {
  if (peptide_idx == 0) {
    return peptide_;
  } else {
    return NULL;
  }
}

/**
 *\returns the number of missed cleavages
 */
int LinearPeptide::getNumMissedCleavages() {
  set<int> skip;
  return peptide_->getMissedCleavageSites(skip);
}

/**
 * \returns whether this peptide is modified by a variable mod
 */
bool LinearPeptide::isModified() {

  return peptide_->isModified();
}

bool compareLinearPeptideMass(
  const LinearPeptide& pep1, 
  const LinearPeptide& pep2) {

  return pep1.getMassConst(MONO) < pep2.getMassConst(MONO);

}
bool compareLinearPeptideMassToFLOAT(const LinearPeptide& pep1, FLOAT_T mass) {

  return pep1.getMassConst(MONO) < mass;
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
