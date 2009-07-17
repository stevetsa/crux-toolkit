#include "hhc.h"

#define NUM_ARGUMENTS 7
#define NUM_OPTIONS 0
using namespace std;

LinkedPeptide::LinkedPeptide(char* peptide_A, char* peptide_B, int posA, int posB, FLOAT_T linkermass, int charge) {
  charge_ = charge;
  linker_mass_ = linkermass;
  Peptide pepA = Peptide(peptide_A);
  Peptide pepB = Peptide(peptide_B);
  pepA.add_link(posA, pepB);
  pepB.add_link(posB, pepA);
  peptides_.push_back(pepA);
  peptides_.push_back(pepB);
}

void LinkedPeptide::calculate_mass(bool bion) {
  mass_ = 0.0;   
  for (vector<Peptide>::iterator peptide = peptides_.begin(); peptide != peptides_.end(); ++peptide) {
    mass_ += calc_sequence_mass((char*)peptide->sequence().c_str(), MONO);
   // change this later!!
  if (peptide != peptides_.begin()) 
    mass_ += linker_mass_;
  } 
  if (bion) 
    mass_ = mass_ - MASS_H2O_MONO;
}

FLOAT_T LinkedPeptide::get_mz() {
  if (mz < 1)
    mz = ((mass_ + MASS_H*charge_) / charge_);
  return mz;
}


void print_ions(vector<LinkedPeptide>& ions) {
  for (vector<LinkedPeptide>::iterator ion = ions.begin(); ion != ions.end(); ++ion) {
    cout << ion->get_mz() << "\t" << "100\t" << *ion << endl;
  }
}

int main(int argc, char** argv) {

  char* peptideA = NULL;
  char* peptideB = NULL;
  int posA = 0;
  int posB = 0;
  int linker_mass = 0;
  int charge = 1; 

 parse_arguments_set_req(
	"peptide A",
 	"first peptide sequence",
	(void *) &peptideA,
	STRING_ARG);

 parse_arguments_set_req(
	"peptide B",
 	"second peptide sequence",
	(void *) &peptideB,
	STRING_ARG);

 parse_arguments_set_req(
	"position A",
 	"zero based index of linker on peptide A [0, length(A)-1]",
	(void *) &posA,
	INT_ARG);

 parse_arguments_set_req(
	"position B",
 	"zero based index of linker on peptide B [0, length(B)-1]",
	(void *) &posB,
	INT_ARG);

 parse_arguments_set_req(
	"charge",
 	"linked peptide charge",
	(void *) &charge,
	INT_ARG);

 parse_arguments_set_req(
	"linker-mass",
 	"mass of the link between A and B",
	(void *) &linker_mass,
	INT_ARG);

  //initialize_parameters(); 
  set_verbosity_level(CARP_INFO);
  //parse_cmd_line_into_params_hash(argc, argv, "");
  //select_cmd_line_arguments( argument_list, num_arguments);
  //select_cmd_line_options( option_list, NUM_OPTIONS );
  parse_arguments(argc, argv, 0);
  carp(CARP_INFO, "(%s:%s),%d,%d, +%d, linker-mass=%d", peptideA, peptideB, posA, posB, charge, linker_mass);  
  LinkedPeptide linked_peptide = LinkedPeptide( peptideA, peptideB, posA, posB, linker_mass, charge);  
  vector<pair<LinkedPeptide, LinkedPeptide> > fragments;
  //linked_peptide.split_many(fragments);
  linked_peptide.split(fragments);
  //print
  vector<LinkedPeptide> all_fragments;

  for (vector<pair<LinkedPeptide, LinkedPeptide> >::iterator it = fragments.begin(); it != fragments.end(); ++it) {
    if (it->first.charge() != 0) {
      it->first.calculate_mass(true); // calculate mass of the b ion
      all_fragments.push_back(it->first);
    }
    if (it->second.charge() != 0) {
      it->second.calculate_mass(false); // calculate mass of y ion
      all_fragments.push_back(it->second);
    }
  }
  for (vector<LinkedPeptide>::iterator lp = all_fragments.begin(); lp != all_fragments.end(); ++lp) {
    lp->get_mz();
  }
  sort(all_fragments.begin(), all_fragments.end());
  print_ions(all_fragments);
  linked_peptide.calculate_mass(false);
  cout << linked_peptide.get_mz() << endl;
  cout << linked_peptide.mass() << endl;
  return 0;
}
