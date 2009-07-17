#include "hhc.h"

#define NUM_ARGUMENTS 7
#define NUM_OPTIONS 0
using namespace std;


// constructor for a single link between two peptides
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

// calculates mass of linked peptide,
// remove H2O from mass if it's a b-ion
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

// 
FLOAT_T LinkedPeptide::get_mz() {
  if (mz < 1)
    mz = ((mass_ + MASS_H*charge_) / charge_);
  return mz;
}

// split the linked peptide at every position, append b and y ions
// to list of all ions
void add_linked_ions(vector<LinkedPeptide>& ions, LinkedPeptide& linked_peptide) {
  vector<pair<LinkedPeptide, LinkedPeptide> > fragments;
  //linked_peptide.split_many(fragments);
  linked_peptide.split(fragments);
  for (vector<pair<LinkedPeptide, LinkedPeptide> >::iterator it = fragments.begin(); it != fragments.end(); ++it) {
    if (it->first.charge() != 0) {
      it->first.calculate_mass(true); // calculate mass of the b ion
      //it->first.get_mz();
      ions.push_back(it->first);
    }
    if (it->second.charge() != 0) {
      it->second.calculate_mass(false); // calculate mass of y ion
      //it->first.get_mz();
      ions.push_back(it->second);
    }
  }
}


void print_ions(vector<LinkedPeptide>& ions) {
  for (vector<LinkedPeptide>::iterator ion = ions.begin(); ion != ions.end(); ++ion) {
    cout << ion->get_mz() << "\t" << "100\t" << *ion << endl;
  }
}

int main(int argc, char** argv) {

  char* fasta_file = NULL;
  char* bonds; 
  int linker_mass = 0;
  int max_charge = 1; 

 parse_arguments_set_req(
	"fasta_file",
 	"database to read peptides from",
	(void *) &fasta_file,
	STRING_ARG);

 parse_arguments_set_req(
	"bonds", 
	"comma delimited list of possible aa linker sites, for example L:A,L:G",
	(void *) &bonds,
	STRING_ARG);

 parse_arguments_set_req(
	"charge",
 	"max possible charge",
	(void *) &max_charge,
	INT_ARG);

 parse_arguments_set_req(
	"linker-mass",
 	"mass of the linker",
	(void *) &linker_mass,
	INT_ARG);

  initialize_parameters(); 
  set_verbosity_level(CARP_INFO);

  parse_arguments(argc, argv, 0);

  // a map to keep track of possible bonds
  string bonds_string = string(bonds);
  map<char, set<char> > bond_map;
  for (int i = 0; i < bonds_string.length() - 2; i += 4) {
     bond_map[bonds_string[i]].insert(bonds_string[i+2]);
  }

  DATABASE_T* database = new_database(fasta_file, FALSE);
  PEPTIDE_CONSTRAINT_T* peptide_constraint = new_peptide_constraint_from_parameters();
  DATABASE_PEPTIDE_ITERATOR_T* peptide_iterator = 
	new_database_peptide_iterator(database, peptide_constraint);
  
 // get a list of pointers to all peptides, using peptide
 // constraint from parameters
 
  vector<PEPTIDE_T*> peptides;
  while (database_peptide_iterator_has_next(peptide_iterator)) {
    peptides.push_back(database_peptide_iterator_next(peptide_iterator));
  }

  // a list of all the ions 
  vector<LinkedPeptide> all_ions;


  // for each peptide in database
  for (vector<PEPTIDE_T*>::iterator peptide = peptides.begin(); peptide != peptides.end(); ++peptide) {
    string alpha_sequence = get_peptide_sequence(*peptide);
    cout << alpha_sequence << endl;
    Peptide alpha = Peptide(alpha_sequence);
    // do something with the single peptide here
    
    // for every other peptide in database
    for (vector<PEPTIDE_T*>::iterator it = peptide+1; it != peptides.end(); ++it) {
      string beta_sequence = get_peptide_sequence(*it);
      Peptide beta = Peptide(beta_sequence);
      // for every aa in alpha
      for (int i = 0; i < alpha_sequence.length(); ++i) {
        map<char, set<char> >::iterator char_it = bond_map.find(alpha_sequence[i]);
        if (char_it != bond_map.end()) { 
	  for (int j = 0; j < beta_sequence.length(); ++j) {
            if (char_it->second.find(beta_sequence[j]) != char_it->second.end()) {
              cout << alpha_sequence << " " << i << " " << beta_sequence << " " << j << endl;
              LinkedPeptide lp = LinkedPeptide((char*) alpha_sequence.c_str(),(char*) beta_sequence.c_str(),i,j,linker_mass,max_charge);
	      add_linked_ions(all_ions, lp);
            } 
	  }
        }
      }
      for (int i = 0; i < beta_sequence.length(); ++i) {
        map<char, set<char> >::iterator char_it = bond_map.find(beta_sequence[i]);
	if (char_it != bond_map.end()) {
	  for (int j = 0; j < alpha_sequence.length(); ++j) {
	    if (char_it->second.find(beta_sequence[j]) != char_it->second.end()) {
              //cout << i << " " << j << endl;
	      LinkedPeptide lp = LinkedPeptide((char*) beta_sequence.c_str(),(char*) alpha_sequence.c_str(),i,j,linker_mass,max_charge);
	      add_linked_ions(all_ions, lp);
	    }
          }
        }
      }
    }
      //cout << all_ions.size() << endl;
  }  
  
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin(); ion != all_ions.end(); ++ion) {
    //ion->calculate_mass(false);
    ion->get_mz();
  }
  sort(all_ions.begin(), all_ions.end());
  print_ions(all_ions);
  return 0;
}
