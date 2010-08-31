#include "DelimitedFile.h"
#include "XLinkBondMap.h"

using namespace std;


const static char* amino_alpha="ABCDEFGHIKLMNPQRSTUVWYZX";
const static int num_amino_alpha = 24;



XLinkBondMap::XLinkBondMap() {
  
  string link_sites = string(get_string_parameter("link sites"));
  setLinkString(link_sites);

}

XLinkBondMap::XLinkBondMap(string& link_sites) {
  setLinkString(link_sites);
}

XLinkBondMap::~XLinkBondMap() {
}

void XLinkBondMap::setLinkString(string& link_sites) {

    //get each bond description
  vector<string> bonds;
  DelimitedFile::tokenize(link_sites, bonds, ',');

  //parse each bond description.
  for (unsigned int bond_idx = 0; bond_idx < bonds.size(); bond_idx++) {
    vector<string> residues;
    DelimitedFile::tokenize(bonds[bond_idx], residues, ':');
    //check for *.
    if (residues[0] == "*" && residues[1] == "*") {
      //add all possible links between residues.
      carp(CARP_INFO, "Linking all residues to each other");
      for (int i=0;i<num_amino_alpha;i++) {
        for (int j=0;j<num_amino_alpha;j++) {
          (*this)[amino_alpha[i]].insert(amino_alpha[j]);
        }
      }
    } else if (residues[0] == "*" || residues[1] == "*") {
      //only one star detected.
      char amino = residues[0][0] == '*' ? residues[1][0] : residues[0][0]; 
      carp(CARP_INFO, "Linking %c to all other residues", amino);
      for (int i=0;i<num_amino_alpha;i++) {
        (*this)[amino].insert(amino_alpha[i]);
        (*this)[amino_alpha[i]].insert(amino);
      }  
    } else {
      //there is no star, insert the link normally
      (*this)[residues[0][0]].insert(residues[1][0]);
      (*this)[residues[1][0]].insert(residues[0][0]);
    }
  }

}

BOOLEAN_T XLinkBondMap::canLink(char aa) {
  return (find(aa) != end());
}

BOOLEAN_T XLinkBondMap::canLink(
  XLinkablePeptide& pep1,
  XLinkablePeptide& pep2,
  int link1_site,
  int link2_site
) {

  char aa1 = get_peptide_sequence_pointer(pep1.getPeptide())[link1_site];
 

  //see if the link exists, this one should always be true.
  XLinkBondMap::iterator find1_iter = find(aa1);
  if (find1_iter == end()) {
    return FALSE;
  }
  //check to see if the second one has a valid link from 1 <-> 2.
  char aa2 = get_peptide_sequence_pointer(pep2.getPeptide())[link2_site];
  set<char>::iterator find2_iter = find1_iter -> second.find(aa2);
  if (find2_iter == find1_iter -> second.end()) {
    return FALSE;
  }

  //sanity check, make sure that the modification doesn't prevent the link.
  //Since an XLinkablePeptide shouldn't contain these, we don't have to
  //recheck here.
  return TRUE;
}

BOOLEAN_T XLinkBondMap::canLink(
  XLinkablePeptide& pep1,
  int link1_site,
  int link2_site) {

  return canLink(pep1, pep1, link1_site, link2_site);
}

BOOLEAN_T XLinkBondMap::canLinkIdx(
  XLinkablePeptide& pep1,
  XLinkablePeptide& pep2,
  int link1_site_idx,
  int link2_site_idx) {
  
  return canLink(pep1, 
		 pep2, 
		 pep1.getLinkSite(link1_site_idx), 
		 pep2.getLinkSite(link2_site_idx));
}
