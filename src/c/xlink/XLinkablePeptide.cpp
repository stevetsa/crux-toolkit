#include "XLinkablePeptide.h"

#include "objects.h"
#include "modifications.h"

#include <iostream>

#include "XLink.h"

using namespace std;

XLinkablePeptide::XLinkablePeptide() {
  peptide_ = NULL;
  sequence_ = NULL;
}

XLinkablePeptide::XLinkablePeptide(char* sequence) {
  cout<<"Assign sequence"<<endl;
  sequence_ = sequence;
  cout<<"Assign peptide"<<endl;
  peptide_ = NULL;
}

XLinkablePeptide::XLinkablePeptide(PEPTIDE_T* peptide,
				   vector<int>& link_sites) {

  peptide_ = peptide;
  link_sites_.clear();
  for (unsigned int idx=0;idx<link_sites.size();idx++) {
    addLinkSite(link_sites[idx]);
  }
}

XLinkablePeptide::XLinkablePeptide(PEPTIDE_T* peptide, 
				   XLinkBondMap& bondmap) {
  peptide_=peptide;
  findLinkSites(peptide_, bondmap, link_sites_);
}

void XLinkablePeptide::findLinkSites(
  PEPTIDE_T* peptide,
  XLinkBondMap& bondmap,
  vector<int>& link_sites) {

  link_sites.clear();

  //find all modifications that can prevent xlink.
  AA_MOD_T** mod_list = NULL;
  int total_mods = get_all_aa_mod_list(&mod_list);
  vector<AA_MOD_T*> prevents_xlink;

  for (int mod_idx = 0;mod_idx < total_mods; mod_idx++) {
    if (aa_mod_get_prevents_xlink(mod_list[mod_idx])) {
      prevents_xlink.push_back(mod_list[mod_idx]);
   }
  }

  //scan through the sequence and find linkable sites.
  MODIFIED_AA_T* mod_seq = get_peptide_modified_aa_sequence(peptide);
  int length = get_peptide_length(peptide);

  char* xlink_prevents_cleavage = 
    get_string_parameter("xlink-prevents-cleavage");
  



  for (int seq_idx=0;seq_idx < length;seq_idx++) {
    if (bondmap.canLink(peptide, seq_idx)) {
      char aa = modified_aa_to_char(mod_seq[seq_idx]);
      BOOLEAN_T link_prevented = FALSE;
      if (seq_idx == length-1) {  
	int idx = 0;
	while (xlink_prevents_cleavage[idx] != '\0') {
	  if (xlink_prevents_cleavage[idx] == aa) {
	    link_prevented = TRUE;
	    break;
	  }
	  idx++;
	}
      }    

      if (!link_prevented) {
	//check if the modification prevents xlink.
	for (unsigned int mod_idx=0;mod_idx<prevents_xlink.size();mod_idx++) {
	  //if aa is modified by any modification that can prevent
	  //cross-linking, then don't add the site as a link_site.
	  if (is_aa_modified(mod_seq[seq_idx], prevents_xlink[mod_idx])) {
	    link_prevented = TRUE;
	    break;
	  }
	}
      }
      //passes both tests, this is a linkable site.
      if (!link_prevented) {
	//if it is a linkable site, then add it to the list.
        link_sites.push_back(seq_idx);
      }
    }
  }
  free(xlink_prevents_cleavage);
  free(mod_seq);
}

XLinkablePeptide::~XLinkablePeptide() {
  //cerr <<"XLinkablePeptide::~XLinkablePeptide()"<<endl;
  //cerr<<"XLinkablePeptide::~XLinkablePeptide():start."<<endl;
  if (peptide_ != NULL) {
    //free_peptide(peptide_);
  }
  //cerr<<"XLinkablePeptide::~XLinkablePeptide():done."<<endl;
}

size_t XLinkablePeptide::numLinkSites() {
  return link_sites_.size();
}

BOOLEAN_T XLinkablePeptide::isLinkable() {
  return numLinkSites() > 0;
}

int XLinkablePeptide::getLinkSite(int link_site_idx) {
  return link_sites_.at(link_site_idx);
}

int XLinkablePeptide::addLinkSite(int link_site) {
  link_sites_.push_back(link_site);
  return link_sites_.size()-1;
}

PEPTIDE_T* XLinkablePeptide::getPeptide() {
  return peptide_;
}

FLOAT_T XLinkablePeptide::getMass() const {
  if (peptide_ == NULL) {
    carp(CARP_INFO,"returning mass for %s", sequence_);
    return calc_sequence_mass(sequence_, AVERAGE);
  } else {
    return get_peptide_peptide_mass(peptide_);
  }
}

FLOAT_T XLinkablePeptide::getMass(MASS_TYPE_T mass_type) {
  if (peptide_ == NULL) {
    return calc_sequence_mass(sequence_, mass_type);
  } else {
    return calc_modified_peptide_mass(peptide_, mass_type);
  }
}


char* XLinkablePeptide::getSequence() {
  if (peptide_ == NULL) {
    return my_copy_string(sequence_);
  } else {
    char* seq = get_peptide_sequence(peptide_);
    return seq;
  }
}

MODIFIED_AA_T* XLinkablePeptide::getModifiedSequence() {
  MODIFIED_AA_T* mod_seq = NULL;

  if (peptide_ == NULL) {
    convert_to_mod_aa_seq(sequence_, &mod_seq);
  } else {
    mod_seq = get_peptide_modified_aa_sequence(peptide_);
  }

  return mod_seq;
}

int findLink(vector<int>& link_sites, int link_site) {
  
  for (unsigned int idx=0;idx<link_sites.size();idx++) {
    if (link_sites[idx] == link_site)
      return idx;
  }
  return -1;

}

void doSwitch(vector<int>& link_sites, int siteA, int siteB) {
  int siteA_link_idx  = findLink(link_sites, siteA);
  int siteB_link_idx = findLink(link_sites, siteB);
  
  if (siteA_link_idx != -1 && siteB_link_idx != -1) {
    int temp_idx = link_sites[siteA_link_idx];
    link_sites[siteA_link_idx] = link_sites[siteB_link_idx];
    link_sites[siteB_link_idx] = temp_idx;
  } else if (siteA_link_idx != -1) {
    link_sites[siteA_link_idx] = siteB;
  } else if (siteB_link_idx != -1) {
    link_sites[siteB_link_idx] = siteA;
  }
}

static const int MAX_SHUFFLES = 5;

char* generateShuffledSequence(
  PEPTIDE_T* peptide,
  vector<int>& link_sites
  ) {

  char* sequence = get_peptide_sequence(peptide);

  int length = get_peptide_length(peptide);

  int num_shuffles = 0;
  do {
    
    //Don't move the n-term and c-term amino acids.
    int start_idx = 1;
    int end_idx = length - 2;

    while (start_idx < end_idx) {
      int switch_idx = get_random_number_interval(start_idx, end_idx);
      char temp_char = sequence[start_idx];
      sequence[start_idx] = sequence[switch_idx];
      sequence[switch_idx] = temp_char;
     
      doSwitch(link_sites, start_idx, switch_idx);

      ++start_idx;
    }


  } while (FALSE /*equal_peptides(sequence, peptide)*/ && (num_shuffles < MAX_SHUFFLES));

  return sequence;


}

MODIFIED_AA_T* generateShuffledModSequence(
  PEPTIDE_T* peptide,
  vector<int>& link_sites
  ) {

  MODIFIED_AA_T* sequence = get_peptide_modified_aa_sequence(peptide);
  int length = get_peptide_length(peptide);
  int start_idx = 0;
  int end_idx = length-1;
  int switch_idx = 0;
  MODIFIED_AA_T temp_aa = 0;

  // Do not move the first and last residue, regardless of enzyme
  ++start_idx;
  --end_idx;
  
  // shuffle from left to right, using the Knuth algorithm for shuffling.
  while(start_idx < end_idx){
    switch_idx = get_random_number_interval(start_idx, end_idx);
    temp_aa = sequence[start_idx];
    sequence[start_idx] = sequence[switch_idx];
    sequence[switch_idx] = temp_aa;
    doSwitch(link_sites, start_idx, switch_idx);
    ++start_idx;
  }

  return sequence;

}


XLinkablePeptide XLinkablePeptide::shuffle() {
  //cerr <<"XLinkablePeptide::shuffle():start"<<endl;
  PEPTIDE_T* peptide = copy_peptide(peptide_);
  XLink::addAllocatedPeptide(peptide);
  //cerr <<"Linksites"<<endl;
  vector<int> link_sites;
  for (unsigned int idx=0;idx<numLinkSites();idx++) {
    link_sites.push_back(getLinkSite(idx));
  }

  MODIFIED_AA_T* new_mod_seq = NULL;

  if(get_peptide_is_modified(peptide)) {
    //cerr<<"Calling generateShuffledModSequence"<<endl;
    new_mod_seq = generateShuffledModSequence(peptide, link_sites);
    set_peptide_decoy_modified_seq(peptide, new_mod_seq);
  } else {
    //cerr<<"Calling generateShuffledSequence"<<endl;
    char* new_seq =
      generateShuffledSequence(peptide, link_sites);
    convert_to_mod_aa_seq(new_seq, &new_mod_seq);
    set_peptide_decoy_modified_seq(peptide, new_mod_seq);
    free(new_seq);
  }
  //cerr<<"Creating new linkable peptide"<<endl;
  XLinkablePeptide ans(peptide, link_sites);
  //cerr <<"XLinkablePeptide::shufle():done."<<endl;
  return ans;
} 

/*
string getAASequence() {
  char* seq = get_peptide_sequence(peptide_);
  string ans(seq);
  free(seq);
  return ans;
}
*/

string XLinkablePeptide::getModifiedSequenceString() {
  if (peptide_ == NULL) {
    return string(sequence_);
  } else {
    char* seq = get_peptide_modified_sequence_with_masses(peptide_, FALSE);
    string string_seq(seq);
    free(seq);
    return string_seq;
  }

}

bool XLinkablePeptide::isModified() {
  return get_peptide_is_modified(peptide_);
}


bool compareXLinkablePeptideMass(
  const XLinkablePeptide& xpep1,
  const XLinkablePeptide& xpep2
) {

  return xpep1.getMass() < xpep2.getMass();
}

