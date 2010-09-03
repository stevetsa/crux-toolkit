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
  linked_peptide_.addLinkSite(posA);
  linked_peptide_.addLinkSite(posB);

  link_pos_idx_.push_back(0);
  link_pos_idx_.push_back(1);

}

SelfLoopPeptide::SelfLoopPeptide(
  XLinkablePeptide& peptide,
  int posA,
  int posB) {

  linked_peptide_ = peptide;
  link_pos_idx_.push_back(posA);
  link_pos_idx_.push_back(posB);
  
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
	    new SelfLoopPeptide(pep, link1_idx, link2_idx);
	  //cerr<<new_candidate -> getSequenceString()<<" " <<new_candidate -> getMass();
	  //cerr<<" "<<min_mass<<" "<<max_mass<<endl;

	  candidates.add(new_candidate);
	}
      }
    }
  }
}


int SelfLoopPeptide::getLinkPos(int link_idx) {
  return linked_peptide_.getLinkSite(link_pos_idx_[link_idx]);
}

MATCHCANDIDATE_TYPE_T SelfLoopPeptide::getCandidateType() {
  return SELFLOOP_CANDIDATE;
}

string SelfLoopPeptide::getSequenceString() {

  string seq = linked_peptide_.getModifiedSequenceString();
  ostringstream oss;
  oss << seq << " (" << (getLinkPos(0)+1) << "," << (getLinkPos(1)+1) << ")";
  string svalue = oss.str();

  return svalue;
}

FLOAT_T SelfLoopPeptide::calcMass(MASS_TYPE_T mass_type) {
  return linked_peptide_.getMass(mass_type) + XLinkPeptide::getLinkerMass();
}

MatchCandidate* SelfLoopPeptide::shuffle() {
  SelfLoopPeptide* decoy = new SelfLoopPeptide();

  decoy->linked_peptide_ = linked_peptide_.shuffle();
  decoy->link_pos_idx_.push_back(link_pos_idx_[0]);
  decoy->link_pos_idx_.push_back(link_pos_idx_[1]);

  return (MatchCandidate*)decoy;


}

void SelfLoopPeptide::predictIons(ION_SERIES_T* ion_series, int charge) {
  char* seq = linked_peptide_.getSequence();
  MODIFIED_AA_T* mod_seq = linked_peptide_.getModifiedSequence();
  set_ion_series_charge(ion_series, charge);
  update_ion_series(ion_series, seq, mod_seq);
  predict_ions(ion_series);
  
  free(mod_seq);

  unsigned int first_site = min(getLinkPos(0), getLinkPos(1));
  unsigned int second_site = max(getLinkPos(0), getLinkPos(1));
  unsigned int N = strlen(seq);
  free(seq);

  //iterate through the ions and modify the ones that have the linker 
  //attached.
  vector<ION_T*> to_remove;
  ION_ITERATOR_T* ion_iter = 
    new_ion_iterator(ion_series);

  while(ion_iterator_has_next(ion_iter)) {
    ION_T* ion = ion_iterator_next(ion_iter);

    bool keep_ion = false;
    bool modify_ion = false;
    unsigned int cleavage_idx = get_ion_cleavage_idx(ion);
   
    if (is_forward_ion_type(ion)) {
      if (cleavage_idx > first_site) {
        if (cleavage_idx <= second_site) {
          keep_ion = false;
        } else {
          keep_ion = true;
          modify_ion = true;
        }
      } else {
        keep_ion = true;
        modify_ion = false;
      }
    } else {
      if (cleavage_idx >= (N-second_site)) {
        if (cleavage_idx < (N-first_site)) {
          keep_ion = false;
        } else {
          keep_ion = true;
          modify_ion = true;
        }
      } else {
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
  

  string ion_sequence = get_ion_peptide_sequence(ion);

  unsigned int cleavage_idx = get_ion_cleavage_idx(ion);

  unsigned int first_site  = min(getLinkPos(0), getLinkPos(1));
  unsigned int second_site = max(getLinkPos(0), getLinkPos(1));

  bool is_linked = false;
  if (is_forward_ion_type(ion)) {
    is_linked = (cleavage_idx > first_site);
  } else {
    is_linked = (cleavage_idx >= (ion_sequence.length() - second_site));
  }

  string subseq;

  //cerr<<"creating substring"<<endl;
  if (is_forward_ion_type(ion)) {
    subseq = ion_sequence.substr(0, cleavage_idx);
  } else {
    subseq = ion_sequence.substr(ion_sequence.length() - cleavage_idx, ion_sequence.length());
  }

  if (is_linked) {
    return subseq + string("*");
  } else {
    return subseq;
  }
}

PEPTIDE_T* SelfLoopPeptide::getPeptide(int peptide_idx) {
  if (peptide_idx == 0) {
    return linked_peptide_.getPeptide();
  } else {
    return NULL;
  }
}
