#include <ctype.h>
#include <unistd.h>
#include "hhc_ion_series.h"
#include "generate_peptides_iterator.h"
#include "parse_arguments.h"
#include "objects.h"
using namespace std;

typedef map<char, set<char> > BondMap;

void find_all_precursor_ions(vector<LinkedPeptide>& all_ions, 
	char* links, 
	FLOAT_T linker_mass, 
	char* missed_link_cleavage,
	char* database_file);

int main(int argc, char** argv) {
  char* min_mass_string = NULL;
  char* max_mass_string = NULL;
  char* database = NULL;
  char* links = NULL;
  char* linker_mass_string = NULL;
  char* missed_link_cleavage = "K";
  int num_missed_cleavages = 0;
  //int charge = 0;
  parse_arguments_set_req(
    "protein database", 
    "database containing all proteins", 
    (void *) &database, 
    STRING_ARG);

  parse_arguments_set_req(
    "links", 
    "comma delimited pair of amino acid link sites, ex. A:K,A:D", 
    (void *) &links, 
    STRING_ARG);
/*
  parse_arguments_set_req(
    "max charge", 
    "maximum charge for ions", 
    (void *) &charge, 
    INT_ARG);
*/
  parse_arguments_set_req(
    "linker mass", 
    "combined mass of linker and linker modifications", 
    (void *) &linker_mass_string, 
    STRING_ARG);
 
  parse_arguments_set_opt(
    "min-mass", 
    "", 
    (void *) &min_mass_string, 
    STRING_ARG);

  parse_arguments_set_opt(
    "max-mass", 
    "", 
    (void *) &max_mass_string, 
    STRING_ARG);

  parse_arguments_set_opt(
    "missed-link-cleavage",
    "",
    (void *) &missed_link_cleavage, 
    STRING_ARG);

  parse_arguments_set_opt(
    "num-missed-cleavages", 
    "maximum number of missed cleavages (not including one at link site)", 
    (void *) &num_missed_cleavages, 
    INT_ARG);

  initialize_parameters();
  parse_arguments(argc, argv, 0);
  // something wrong with DOUBLE_ARG
  FLOAT_T linker_mass = atof(linker_mass_string);
  
  //atof(min_mass_string);
  //atof(max_mass_string);

  FLOAT_T c = get_double_parameter("C"); 
  cout << "c " << c << endl;
  vector<LinkedPeptide> all_ions;
  
  find_all_precursor_ions(all_ions, links, linker_mass, missed_link_cleavage, database);
  
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin(); ion != all_ions.end(); ++ion) {
    ion->calculate_mass();
  }
 
  sort(all_ions.begin(), all_ions.end());

  FLOAT_T max_mass = all_ions.back().mass();
  //FLOAT_T max_mass = 1400;
  FLOAT_T min_mass = 0.0;
  if (min_mass_string != NULL) min_mass = atof(min_mass_string);
  if (max_mass_string != NULL) max_mass = atof(max_mass_string);
  if (max_mass < min_mass) {
    carp(CARP_FATAL, "max mass must be larger than min mass");
  }

  cout << "min " << min_mass << " max " << max_mass << endl;

  LinkedIonSeries ion_fragments = LinkedIonSeries();
   
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin(); ion != all_ions.end(); ++ion) {
    if (min_mass <= ion->mass() && ion->mass() <= max_mass) {
      cout << ion->mass() << "\t" << *ion << endl;
      ion_fragments.add_linked_ions(*ion);
    }
  }
  
  vector<LinkedPeptide> fragments = ion_fragments.ions();
  //cout << "linked precursors: " << ion_fragments.size() << endl;

  sort(fragments.begin(), fragments.end());
  
  for (vector<LinkedPeptide>::iterator fragment = fragments.begin();
		fragment != fragments.end(); ++fragment) {
    //cout << fragment->mass() << "\t" << *fragment << endl;
  } 
 
   
  cout << "links: " << links << endl;
  cout << "linker mass: " << linker_mass << endl;
  //cout << "mass\tpeptide" << endl;
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin(); ion != all_ions.end(); ++ion) {
    //cout << ion->mass() << "\t" << *ion << endl;
  }

  cout << "ions: " << all_ions.size() << endl;
   
  return 0;
}

// a hack, works for EDC linker only
void get_linkable_peptides(set<string>& peptides, 
	DATABASE_PROTEIN_ITERATOR_T* protein_iterator,
	PEPTIDE_CONSTRAINT_T* peptide_constraint) 
{
  PROTEIN_PEPTIDE_ITERATOR_T* peptide_iterator = NULL;
  PROTEIN_T* protein;
  PEPTIDE_T* peptide;
  string sequence = "";
  string last_sequence = "zz";
  bool missed_cleavage = false;
  // keep track of whether the next peptide contains the previous one or not
  size_t index;
  while (database_protein_iterator_has_next(protein_iterator)) {
    protein = database_protein_iterator_next(protein_iterator);
    peptide_iterator = new_protein_peptide_iterator(protein, peptide_constraint);
    // missed_cleavages must be TRUE in protein.c for this to work
    prepare_protein_peptide_iterator(peptide_iterator); 
    while (protein_peptide_iterator_has_next(peptide_iterator)) {
      //peptide = database_peptide_iterator_next(peptide_iterator);
      peptide = protein_peptide_iterator_next(peptide_iterator);
      sequence = get_peptide_sequence(peptide); 
      index = sequence.find(last_sequence);
      // if doesn't contain last peptide
      if (sequence[0] == 'R' && sequence[1] != 'P') { continue;}
      if (index == string::npos || missed_cleavage) {
        missed_cleavage = !missed_cleavage;
        if (!missed_cleavage && last_sequence[last_sequence.size()-1] != 'K') {
          //cout << "skipping " << get_peptide_sequence(peptide) << endl;
	  continue;
	} 
        cout << "peptide " << get_peptide_sequence(peptide) << endl;
        peptides.insert(string(get_peptide_sequence(peptide)));
      } else {
	//cout << "skipping " << get_peptide_sequence(peptide) << endl;
	missed_cleavage = false;
      }
      last_sequence = string(get_peptide_sequence(peptide));
    }
  } 
}

void find_all_precursor_ions(vector<LinkedPeptide>& all_ions, 
			     char* links, 
			     FLOAT_T linker_mass,
			     char* missed_link_cleavage,
		             char* database_file)
{
 
  DATABASE_T* db = new_database(database_file, FALSE);
  PEPTIDE_CONSTRAINT_T* peptide_constraint = new_peptide_constraint_from_parameters();
  // add 
  set_peptide_constraint_num_mis_cleavage(peptide_constraint, 1);
  //set_verbosity_level(CARP_INFO);
  PROTEIN_T* protein = NULL;
  DATABASE_PROTEIN_ITERATOR_T* protein_iterator = new_database_protein_iterator(db);
  PROTEIN_PEPTIDE_ITERATOR_T* peptide_iterator = NULL;

  BondMap bonds; 
  string bonds_string = string(links);
  for (int i = 0; i < bonds_string.length() - 2; i += 4) {
     bonds[bonds_string[i+2]].insert(bonds_string[i]);
     bonds[bonds_string[i]].insert(bonds_string[i+2]);
  }
 
  set<string> peptides;
  get_linkable_peptides(peptides, protein_iterator, peptide_constraint);
  cout << "peptides: " << peptides.size() << endl;
   
  for (set<string>::iterator pepA = peptides.begin(); pepA != peptides.end(); ++pepA) {
    char* sequenceA = (char*) pepA->c_str();
    // add unlinked precursor
    LinkedPeptide lp = LinkedPeptide(1, linker_mass);
    Peptide p = Peptide(sequenceA);
    lp.add_peptide(p);
    all_ions.push_back(lp);
    for (int i = 0; i < pepA->length(); ++i) {
      BondMap::iterator bond = bonds.find(pepA->at(i));
      // if a link aa and doesn't end in K
      if (bond != bonds.end() && i != pepA->length()-1) {
	if (i == pepA->length()-1 && pepA->at(pepA->length()-1) == 'K') continue;
        // add dead end
	all_ions.push_back(LinkedPeptide(sequenceA, NULL, i, -1, linker_mass, 1));
        // add self loop
	for (int j = i+1; j < pepA->length(); ++j) {
          if (bond->second.find(pepA->at(j)) != bond->second.end()) { 
	    //skip if linked to a K at the end
	    if (j == pepA->length()-1 && pepA->at(pepA->length()-1) == 'K') continue;
	    all_ions.push_back(LinkedPeptide(sequenceA, NULL, i, j, linker_mass,1));
	  }
	}
        // add linked precursor
        for (set<string>::iterator pepB = pepA; pepB != peptides.end(); ++pepB) {
          char* sequenceB = (char*) pepB->c_str();
          for (int j = 0; j < pepB->length(); ++j) {
	    // if a link aa and doesn't end in K
	    if (bond->second.find(pepB->at(j)) != bond->second.end()) { 
	      // skip if link to K at end
	      if (j == pepB->length()-1 && pepB->at(pepB->length()-1) == 'K') continue;
	      all_ions.push_back(LinkedPeptide(sequenceA, sequenceB, i, j, linker_mass,1));
	    }
          }
	} // get next pepB 
      }
    }
  } // get next pepA 
} 
 
