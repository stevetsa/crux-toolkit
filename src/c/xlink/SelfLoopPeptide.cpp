#include "SelfLoopPeptide.h"
#include "XLinkablePeptide.h"
#include "XLinkPeptide.h"
#include "XLink.h"

#include "IonSeries.h"
#include "Ion.h"

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

void SelfLoopPeptide::addCandidates(
  FLOAT_T min_mass,
  FLOAT_T max_mass,
  XLinkBondMap& bondmap, 
  Index* index, Database* database,
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods,
  XLinkMatchCollection& candidates) {

  //int max_missed_cleavages = get_int_parameter("missed-cleavages");
  
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
	if (bondmap.canLink(pep, link1_idx, link2_idx)) {
	  //create the candidate.
	  //cerr<<"Adding new selfloop peptide"<<endl;
	  XLinkMatch* new_candidate = 
	    new SelfLoopPeptide(pep, link1_idx, link2_idx);

          //if (new_candidate->getNumMissedCleavages() <= max_missed_cleavages) {
            candidates.add(new_candidate);
          //} else {
          //  delete new_candidate;
          //}
	}
      }
    }
  }
}


int SelfLoopPeptide::getLinkPos(int link_idx) {
  return linked_peptide_.getLinkSite(link_pos_idx_[link_idx]);
}

XLINKMATCH_TYPE_T SelfLoopPeptide::getCandidateType() {
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

XLinkMatch* SelfLoopPeptide::shuffle() {
  SelfLoopPeptide* decoy = new SelfLoopPeptide();

  decoy->linked_peptide_ = linked_peptide_.shuffle();
  decoy->link_pos_idx_.push_back(link_pos_idx_[0]);
  decoy->link_pos_idx_.push_back(link_pos_idx_[1]);

  return (XLinkMatch*)decoy;


}

void SelfLoopPeptide::predictIons(IonSeries* ion_series, int charge) {
  char* seq = linked_peptide_.getSequence();
  MODIFIED_AA_T* mod_seq = linked_peptide_.getModifiedSequence();
  ion_series->setCharge(charge);
  ion_series->update(seq, mod_seq);
  ion_series->predictIons();
  
  free(mod_seq);

  unsigned int first_site = min(getLinkPos(0), getLinkPos(1));
  unsigned int second_site = max(getLinkPos(0), getLinkPos(1));
  unsigned int N = strlen(seq);
  free(seq);

  //iterate through the ions and modify the ones that have the linker 
  //attached.
  vector<Ion*> to_remove;

  for (IonIterator ion_iter = ion_series->begin();
    ion_iter != ion_series->end();
    ++ion_iter) {

    Ion* ion = *ion_iter;

    bool keep_ion = false;
    bool modify_ion = false;
    unsigned int cleavage_idx = ion->getCleavageIdx();
   
    if (ion->isForwardType()) {
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
	FLOAT_T mass_z = ion->getMassZ();
	int charge = ion->getCharge();
	double mass = (mass_z -MASS_PROTON) * (double)charge;
	mass += XLinkPeptide::getLinkerMass();
	mass_z = (mass + MASS_PROTON * (double)charge) / (double)charge;
	ion->setMassZ(mass_z);
      }
    } else {
      to_remove.push_back(ion);
    }
  }

  for (unsigned int idx=0;idx < to_remove.size();idx++) {
    ion_series->removeIon(to_remove[idx]);
    delete to_remove[idx];
  }

}

string SelfLoopPeptide::getIonSequence(Ion* ion) {
  

  string ion_sequence = ion->getPeptideSequence();

  unsigned int cleavage_idx = ion->getCleavageIdx();

  unsigned int first_site  = min(getLinkPos(0), getLinkPos(1));
  unsigned int second_site = max(getLinkPos(0), getLinkPos(1));

  bool is_linked = false;
  if (ion->isForwardType()) {
    is_linked = (cleavage_idx > first_site);
  } else {
    is_linked = (cleavage_idx >= (ion_sequence.length() - second_site));
  }

  string subseq;

  //cerr<<"creating substring"<<endl;
  if (ion->isForwardType()) {
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

int SelfLoopPeptide::getNumMissedCleavages() {
  char missed_cleavage_link_site = 'K';

  int link1_site = getLinkPos(0);
  int link2_site = getLinkPos(1);

  set<int> skip;

  PEPTIDE_T* pep = linked_peptide_.getPeptide();

  char* seq = get_peptide_sequence_pointer(pep);

  if (seq[link1_site] == missed_cleavage_link_site) {
    skip.insert(link1_site);
  }

  if (seq[link2_site] == missed_cleavage_link_site) {
    skip.insert(link2_site);
  }

  return get_peptide_missed_cleavage_sites(pep, skip);

}

bool SelfLoopPeptide::isModified() {

  return linked_peptide_.isModified();
}
