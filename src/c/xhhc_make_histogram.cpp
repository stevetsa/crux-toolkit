#include "hhc_ion_series.h"
#include "objects.h"

using namespace std;

typedef map<char, set<char> > BondMap;

void find_all_precursor_ions(vector<LinkedPeptide>& all_ions, 
	char* links, 
	FLOAT_T linker_mass, 
	char* database_file);

int main(int argc, char** argv) {
  char* ms2_file = NULL;
  char* min_mass_string = NULL;
  char* max_mass_string = NULL;
  char* database = NULL;
  char* links = NULL;
  char* linker_mass_string = NULL;
  int charge = 1;
  int scan_num = 0;
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

  parse_arguments_set_req(
    "scan-number", 
    "The scan number for the MS-MS spectrum to extract from the ms2 file. This is an integer in the range [1, 100000], and uniquely identifies a particular MS-MS spectrum within an .ms2 file.",
    (void *) &scan_num, INT_ARG);

  parse_arguments_set_req(
    "ms2-filename", 
    "A file containing multiple MS-MS spectra in .ms2 format.",
    (void *) &ms2_file,
    STRING_ARG);
 
  parse_arguments_set_opt(
    "charge",
    "peptide charge", 
    (void *) &charge,
    INT_ARG); 

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

  initialize_parameters();
  parse_arguments(argc, argv, 0);

  // something wrong with DOUBLE_ARG
  FLOAT_T linker_mass = atof(linker_mass_string);
  
  //atof(min_mass_string);
  //atof(max_mass_string);

  cout << "ms2 " << ms2_file << " charge " << charge << " scan num " << scan_num << endl;

  vector<LinkedPeptide> all_ions;
  
  find_all_precursor_ions(all_ions, links, linker_mass, database);
  
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin(); ion != all_ions.end(); ++ion) {
    ion->calculate_mass();
  }
 
  sort(all_ions.begin(), all_ions.end());

  FLOAT_T max_mass = all_ions.back().mass();
  FLOAT_T min_mass = 0.0;
  if (min_mass_string != NULL) min_mass = atof(min_mass_string);
  if (max_mass_string != NULL) max_mass = atof(max_mass_string);
  if (max_mass < min_mass) {
    carp(CARP_FATAL, "max mass must be larger than min mass");
  }

  cout << "min " << min_mass << " max " << max_mass << endl;

  //LinkedIonSeries ion_fragments = LinkedIonSeries();
  vector<LinkedPeptide> filtered_ions;
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin(); ion != all_ions.end(); ++ion) {
    if (min_mass <= ion->mass() && ion->mass() <= max_mass) {
      ion->set_charge(charge);
      filtered_ions.push_back(*ion);
      cout << ion->mass() << "\t" << *ion << endl;
    } 
      
  }
  
  //vector<LinkedPeptide> fragments = ion_fragments.ions();
  cout << filtered_ions.size() << endl;
  
   
    SPECTRUM_T* spectrum = NULL;
    SPECTRUM_COLLECTION_T* collection = NULL;
    SCORER_T* scorer = NULL;
    FLOAT_T score = 0;
    LinkedIonSeries ion_series;
    map<FLOAT_T, int> buckets;
    int count = 0;
  for (vector<LinkedPeptide>::iterator frag = filtered_ions.begin(); frag != filtered_ions.end(); ++frag) {
    ion_series = LinkedIonSeries(links, charge, linker_mass);
    ion_series.add_linked_ions(*frag);
    collection = new_spectrum_collection(ms2_file);
    spectrum = allocate_spectrum();
    scorer = new_scorer(XCORR);
    vector<LinkedPeptide> series = ion_series.ions();

    if(!get_spectrum_collection_spectrum(collection, scan_num, spectrum)){
      carp(CARP_ERROR, "failed to find spectrum with  scan_num: %d", scan_num);
      free_spectrum_collection(collection);
      free_spectrum(spectrum);
      exit(1);
    }

    score = hhc_score_spectrum_v_ion_series(scorer, spectrum, ion_series);
    //cout << count << " " << floor(score*10)/10 << endl;
    ++buckets[floor(score*10) / 10]; 
    free_scorer(scorer);
    free_spectrum_collection(collection);
    free_spectrum(spectrum);
    ++count;
  }
/*
  sort(fragments.begin(), fragments.end());
  
  for (vector<LinkedPeptide>::iterator fragment = fragments.begin();
		fragment != fragments.end(); ++fragment) {
    cout << fragment->mass() << "\t" << *fragment << endl;
  } 
 */
  /* 
  cout << "links: " << links << endl;
  cout << "linker mass: " << linker_mass << endl << endl;
  cout << "mass\tpeptide" << endl;
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin(); ion != all_ions.end(); ++ion) {
    cout << ion->mass() << "\t" << *ion << endl;
  }
*/
  cout << "ions: " << all_ions.size() << endl;
  for (map<FLOAT_T, int>::iterator it = buckets.begin(); it != buckets.end(); ++it) {
    cout << it->first << " " << it->second << endl;
  }  

  return 0;
}

void find_all_precursor_ions(vector<LinkedPeptide>& all_ions, 
			     char* links, 
			     FLOAT_T linker_mass,
		             char* database_file)
{
  DATABASE_T* db = new_database(database_file, FALSE);
  PEPTIDE_CONSTRAINT_T* peptide_constraint = new_peptide_constraint_from_parameters();
  //set_verbosity_level(CARP_INFO);

  DATABASE_PEPTIDE_ITERATOR_T* peptide_iterator = 
	new_database_peptide_iterator(db, peptide_constraint);
  BondMap bonds; 
  string bonds_string = string(links);
  for (int i = 0; i < bonds_string.length() - 2; i += 4) {
     bonds[bonds_string[i+2]].insert(bonds_string[i]);
     bonds[bonds_string[i]].insert(bonds_string[i+2]);
  }
 
  vector<string> peptides;
  PEPTIDE_T* peptide = NULL; 
  while (database_peptide_iterator_has_next(peptide_iterator)) {
    peptide = database_peptide_iterator_next(peptide_iterator);
    peptides.push_back(string(get_peptide_sequence(peptide)));
  } 
 
  for (vector<string>::iterator pepA = peptides.begin(); pepA != peptides.end(); ++pepA) {
    char* sequenceA = (char*) pepA->c_str();
    // add unlinked precursor
    LinkedPeptide lp = LinkedPeptide(1, linker_mass);
    Peptide p = Peptide(sequenceA);
    lp.add_peptide(p);
    all_ions.push_back(lp);
    for (int i = 0; i < pepA->length(); ++i) {
      BondMap::iterator bond = bonds.find(pepA->at(i));
      if (bond != bonds.end()) {
        // add dead end
	all_ions.push_back(LinkedPeptide(sequenceA, NULL, i, -1, linker_mass, 1));
        // add self loop
	for (int j = i+1; j < pepA->length(); ++j) {
          if (bond->second.find(pepA->at(j)) != bond->second.end()) 
	    all_ions.push_back(LinkedPeptide(sequenceA, NULL, i, j, linker_mass,1));
	}
        // add linked precursor
        for (vector<string>::iterator pepB = pepA + 1; pepB != peptides.end(); ++pepB) {
          char* sequenceB = (char*) pepB->c_str();
          for (int j = 0; j < pepB->length(); ++j) {
	    if (bond->second.find(pepB->at(j)) != bond->second.end()) 
	      all_ions.push_back(LinkedPeptide(sequenceA, sequenceB, i, j, linker_mass,1));
          }
	} // get next pepB 
      }
    }
  } // get next pepA
} 
 
