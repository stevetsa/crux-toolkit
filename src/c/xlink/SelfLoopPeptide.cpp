#include "SelfLoopPeptide.h"
#include "XLinkablePeptide.h"
#include "XLinkPeptide.h"
#include "XLink.h"

#include "ion_series.h"
#include "ion.h"

#include <iostream>
#include <sstream>

using namespace std;

SelfLoopPeptide::SelfLoopPeptide() {
}

SelfLoopPeptide::SelfLoopPeptide(
  char* peptide,
  int posA,
  int posB) {

  linked_peptide_ = XLinkablePeptide(peptide);
  link_pos_.push_back(posA);
  link_pos_.push_back(posB);

  sort(link_pos_.begin(), link_pos_.end());

}

SelfLoopPeptide::SelfLoopPeptide(
  XLinkablePeptide& peptide,
  int posA,
  int posB) {

  linked_peptide_ = peptide;
  link_pos_.push_back(posA);
  link_pos_.push_back(posB);

  sort(link_pos_.begin(), link_pos_.end());
  
}

SelfLoopPeptide::~SelfLoopPeptide() {
}

void SelfLoopPeptide::addCandidates(FLOAT_T precursor_mz, int charge, 
  XLinkBondMap& bondmap, 
  INDEX_T* index, DATABASE_T* database,
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods,
  MatchCandidateVector& candidates,
  BOOLEAN_T use_decoy_window) {

  FLOAT_T min_mass,max_mass;

  get_min_max_mass(precursor_mz, 
		   charge,
		   use_decoy_window,
		   min_mass,
		   max_mass);
  
  vector<XLinkablePeptide> linkable_peptides;
  int cur_aa_mods = 0;
  //loop over modifications.
  for (int mod_idx=0;mod_idx<num_peptide_mods;mod_idx++) {
  
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];
    int this_aa_mods = peptide_mod_get_num_aa_mods(peptide_mod);
    if (this_aa_mods > cur_aa_mods) {
      cur_aa_mods = this_aa_mods;
    }
    XLinkPeptide::addLinkablePeptides(
      min_mass - XLinkPeptide::getLinkerMass(),
      max_mass - XLinkPeptide::getLinkerMass(),
      index,
      database,
      peptide_mod,
      FALSE,
      bondmap,
      linkable_peptides);
  }
  //cerr<<"We have "<<linkable_peptides.size()<<" linkable peptides"<<endl;
  //find linkable peptides that can have links to themselves.
  for (unsigned int idx =0;idx < linkable_peptides.size();idx++) {
    XLinkablePeptide &pep = linkable_peptides[idx];

    for (unsigned int link1_idx=0;link1_idx<pep.numLinkSites()-1;link1_idx++) {
      for (unsigned int link2_idx=link1_idx+1;link2_idx<pep.numLinkSites();link2_idx++) {
	int link1_site = pep.getLinkSite(link1_idx);
	int link2_site = pep.getLinkSite(link2_idx);
	if (bondmap.canLink(pep, link1_site, link2_site)) {
	  //create the candidate.
	  //cerr<<"Adding new selfloop peptide"<<endl;
	  MatchCandidate* new_candidate = 
	    new SelfLoopPeptide(pep, link1_site, link2_site);
	  //cerr<<new_candidate -> getSequenceString()<<" " <<new_candidate -> getMass();
	  //cerr<<" "<<min_mass<<" "<<max_mass<<endl;

	  candidates.add(new_candidate);
	}
      }
    }
  }
}

MATCHCANDIDATE_TYPE_T SelfLoopPeptide::getCandidateType() {
  return SELFLOOP_CANDIDATE;
}

string SelfLoopPeptide::getSequenceString() {
  char* seq = get_peptide_modified_sequence_with_masses(linked_peptide_.getPeptide(), FALSE);
  
  ostringstream oss;

  oss << seq << " (" << (link_pos_[0]+1) << "," << (link_pos_[1]+1) << ")";

  string svalue = oss.str();

  free(seq);

  return svalue;
}

FLOAT_T SelfLoopPeptide::getMass() {
  return linked_peptide_.getMass() + XLinkPeptide::getLinkerMass();
}

MatchCandidate* SelfLoopPeptide::shuffle() {
  SelfLoopPeptide* decoy = new SelfLoopPeptide(*this);


  return (MatchCandidate*)decoy;


}

void SelfLoopPeptide::predictIons(ION_SERIES_T* ion_series, int charge) {
  char* seq = get_peptide_sequence(linked_peptide_.getPeptide());
  MODIFIED_AA_T* mod_seq = get_peptide_modified_aa_sequence(linked_peptide_.getPeptide());
  set_ion_series_charge(ion_series, charge);
  update_ion_series(ion_series, seq, mod_seq);
  predict_ions(ion_series);
  free(seq);
  free(mod_seq);


  //iterate through the ions and modify the ones that have the linker 
  //attached.
  vector<ION_T*> to_remove;
  ION_ITERATOR_T* ion_iter = 
    new_ion_iterator(ion_series);

  while(ion_iterator_has_next(ion_iter)) {
    int start_idx=0;
    int end_idx=0;
    ION_T* ion = ion_iterator_next(ion_iter);

    bool keep_ion = false;
    bool modify_ion = false;

    switch(get_ion_type(ion)) {
    case B_ION:
      start_idx = 0;
      end_idx = get_ion_cleavage_idx(ion)+1; 
      if (link_pos_[0] == 0) {
	if (end_idx >= link_pos_[1]) {
	  keep_ion = true;
	  modify_ion = true;
	}
      }
      break;
    case Y_ION:
      start_idx = get_ion_cleavage_idx(ion);
      end_idx = get_peptide_length(linked_peptide_.getPeptide());
      if (link_pos_[1] == end_idx) {
	if (start_idx <= link_pos_[0]) {
	  keep_ion = true;
	  modify_ion = true;
	}
      }
    default:
      keep_ion = false;
      modify_ion = false;
    }

    if (!keep_ion) {
      if (start_idx < link_pos_[0] && end_idx < link_pos_[0]) {
	keep_ion = true;
	modify_ion = false;
      } else if (start_idx > link_pos_[1] && end_idx > link_pos_[1]) {
	keep_ion = true;
	modify_ion = false;
      }
    }

    if (keep_ion) {
      if (modify_ion) {
	FLOAT_T mass_z = get_ion_mass_z(ion);
	int charge = get_ion_charge(ion);
	double mass = (mass_z -MASS_PROTON) * (double)charge;
	mass += XLinkPeptide::getLinkerMass();
	mass_z = (mass + MASS_PROTON * (double)charge) / (double)charge;
	set_ion_mass_z(ion, mass_z);
      }
    } else {
      to_remove.push_back(ion);
    }


  }

  for (unsigned int idx=0;idx < to_remove.size();idx++) {
    remove_ion_from_ion_series(ion_series, to_remove[idx]);
    free_ion(to_remove[idx]);
  }

  free_ion_iterator(ion_iter);



}

string SelfLoopPeptide::getIonSequence(ION_T* ion) {
  return get_ion_peptide_sequence(ion);
}
