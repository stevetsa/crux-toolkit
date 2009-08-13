#include "hhc_ion_series.h"
#include "objects.h"
#include <fstream>
#include <math.h>
//#include <iostream>

//#define PARAM_ESTIMATION_SAMPLE_COUNT 500
#define MIN_WEIBULL_MATCHES 40
#define MIN_XCORR_SHIFT -5.0
#define MAX_XCORR_SHIFT  5.0
#define XCORR_SHIFT 0.05
using namespace std;

typedef map<char, set<char> > BondMap;

BOOLEAN_T hhc_estimate_weibull_parameters_from_xcorrs(
  FLOAT_T* scores,
  int num_scores,
  FLOAT_T* eta,
  FLOAT_T* beta,
  FLOAT_T* shift,
  FLOAT_T* correlation,
  SPECTRUM_T* spectrum,
  int charge
  );

void add_linked_peptides(vector<LinkedPeptide>& all_ions, set<string>& peptides, string links, FLOAT_T linker_mass, int charge, bool is_decoy);

//void add_decoys(vector<LinkedPeptide>& decoys, LinkedPeptide& lp, char* links, int charge, FLOAT_T linker_mass);

void add_decoys(vector<LinkedPeptide>& decoys, LinkedPeptide& lp);

void plot_weibull(vector<pair<FLOAT_T, LinkedPeptide> >& scores, SPECTRUM_T* spectrum, int charge); 

void find_all_precursor_ions(vector<LinkedPeptide>& all_ions, 
	char* links, 
	FLOAT_T linker_mass,
	int charge,
	char* missed_link_cleavage, 
	char* database_file);

int main(int argc, char** argv) {
  char* missed_link_cleavage = "K";
  int num_missed_cleavages = 0;
  char* ms2_file = NULL;
  char* min_mass_string = NULL;
  char* max_mass_string = NULL;
  char* database = NULL;
  char* links = NULL;
  char* linker_mass_string = NULL;
  int decoy_iterations = 5;
  int charge = 1;
  int scan_num = 0;
  bool open_modification = false;
  int open_modification_int = 0;
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
    "open-modification",
    "",
    (void *) &open_modification_int, 
    INT_ARG);

  parse_arguments_set_opt(
    "decoy-iterations",
    "",
    (void *) &decoy_iterations,
    INT_ARG);

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
  // and boolean arg
  if (open_modification_int == 1)
    open_modification = true;
  cout << "ms2 " << ms2_file << " charge " << charge << " scan num " << scan_num << endl;

  vector<LinkedPeptide> all_ions;
  
  find_all_precursor_ions(all_ions, links, linker_mass, charge, missed_link_cleavage, database);

  FLOAT_T max_mass = all_ions.back().mass();
  FLOAT_T min_mass = 0.0;
  if (min_mass_string != NULL) min_mass = atof(min_mass_string);
  if (max_mass_string != NULL) max_mass = atof(max_mass_string);
  if (max_mass < min_mass) {
    carp(CARP_FATAL, "max mass must be larger than min mass");
  }

  cout << "min " << min_mass << " max " << max_mass << endl;

  int num_ions = 0;
  vector<LinkedPeptide> filtered_ions;
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin(); ion != all_ions.end(); ++ion) {
      ion->calculate_mass();
    // if the mass is in the range
    if (min_mass <= ion->mass() && ion->mass() <= max_mass) {
      //ion->set_charge(charge);
      ++num_ions;
      filtered_ions.push_back(*ion);
      // print out the ion
      cout << ion->mass() << "\t" << *ion << endl;
      // iterate to add shuffled decoys
      //int i = 5;
      //if (charge == 4) i = 7;
      for (int i = decoy_iterations; i > 0; --i)
        add_decoys(filtered_ions, *ion);
        //add_decoys(filtered_ions, *ion, links, charge, linker_mass);
    } 
 }     
  
  // sort filtered ions and decoy ions by mass
  cout << "sorting ...";
  sort(filtered_ions.begin(), filtered_ions.end());
  cout << "done" << endl;
  cout << "scan " << scan_num << " +" << charge << "<br>" << endl;
  cout << "precursors  " << num_ions << "<br>" << endl;
  cout << "decoys      " << filtered_ions.size() - num_ions << "<br>" << endl;
  cout << "total       " << filtered_ions.size() << "<br>" << endl;

  SPECTRUM_T* spectrum = allocate_spectrum();
  SPECTRUM_COLLECTION_T* collection = new_spectrum_collection(ms2_file);
  SCORER_T* scorer = new_scorer(XCORR);
  
  if(!get_spectrum_collection_spectrum(collection, scan_num, spectrum)){
    carp(CARP_ERROR, "failed to find spectrum with  scan_num: %d", scan_num);
    free_spectrum_collection(collection);
    free_spectrum(spectrum);
    exit(1);
  }
  
  FLOAT_T score = 0;
  FLOAT_T mod_mass;
  // Pragya's open modification method
  if (open_modification) {
    ION_SERIES_T* ion_series = NULL; 
    ION_CONSTRAINT_T* ion_constraint = 
	new_ion_constraint_sequest_xcorr(charge);
    set<pair<FLOAT_T, string> > scores;
    stringstream ss;
    // for every precursor in the mass window
    for (vector<LinkedPeptide>::iterator ion = filtered_ions.begin(); ion != filtered_ions.end(); ++ion) {
      if (ion->size() == 2) {
	vector<Peptide> peptides = ion->peptides();
	// score the first peptide with modification of second peptide	
        mod_mass = linker_mass + peptides[1].mass();	
	ion_series = new_ion_series((char*)peptides[0].sequence().c_str(), ion->charge(), ion_constraint);
	hhc_predict_ions(ion_series, mod_mass, peptides[0].link_site());
	score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);
	ss.str("");
	ss << peptides[0].sequence() << " mod " << peptides[1].sequence() << ", " << peptides[0].link_site();
        scores.insert(make_pair(score, ss.str()));	
	// score second peptide with modification of first peptide
        mod_mass = linker_mass + peptides[0].mass();	
	ion_series = new_ion_series((char*)peptides[1].sequence().c_str(), ion->charge(), ion_constraint);
	hhc_predict_ions(ion_series, mod_mass, peptides[1].link_site());
	score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);
	ss.str("");
	ss << peptides[1].sequence() << " mod " << peptides[0].sequence() << ", " << peptides[1].link_site();
        scores.insert(make_pair(score, ss.str()));
      }
    }
    //sort(scores.begin(), scores.end());
    int i = 1;
    for (set<pair<FLOAT_T, string> >::reverse_iterator score_pair = scores.rbegin();
		score_pair != scores.rend();
		++score_pair) {

	cout << i << "\t" << score_pair->first << "\t" << score_pair->second << endl;
	++i;
    }

  } else { // linked peptide method
    vector<pair<FLOAT_T, LinkedPeptide> > scores;
    //LinkedIonSeries ion_series;
  // for every ion in the mass window
    for (vector<LinkedPeptide>::iterator ion = filtered_ions.begin(); ion != filtered_ions.end(); ++ion) {
      LinkedIonSeries ion_series = LinkedIonSeries(links, charge, linker_mass);
      ion_series.add_linked_ions(*ion);
      vector<LinkedPeptide> series = ion_series.ions();
      score = hhc_score_spectrum_v_ion_series(scorer, spectrum, ion_series);
      scores.push_back(make_pair(score, *ion));
    }
    sort(scores.begin(), scores.end());
    FLOAT_T range = scores.back().first - scores.front().first;
    cout << "xcorr range " << range << "<br>" << endl;
    plot_weibull(scores, spectrum, charge);
    reverse(scores.begin(), scores.end());
    for (int i = 0; i < 20; ++i) {
	cout << scores[i].first << "\t" << scores[i].second << endl;
    }
  }
  free_scorer(scorer);
  free_spectrum_collection(collection);
  free_spectrum(spectrum);
  return 0;
}

// for running experiments. plots fit and pvalues
void plot_weibull(vector<pair<FLOAT_T, LinkedPeptide> >& scores, SPECTRUM_T* spectrum, int charge) {
  
  ofstream target_fit_file ("fit.target");
  ofstream decoy_fit_file ("fit.decoy");
  ofstream target_score_file ("scores.target");
  ofstream decoy_score_file ("scores.decoy");
  ofstream target_pvalue_file ("pvalues.target");
  ofstream decoy_pvalue_file ("pvalues.decoy");

  int num_scores = scores.size();
  int num_targets = 0;
  int num_decoys = 0;

  FLOAT_T decoy_scores_array[num_scores];
  FLOAT_T target_scores_array[num_scores];

  for (vector<pair<FLOAT_T, LinkedPeptide> >::iterator score_pair = scores.begin();
	score_pair != scores.end(); ++score_pair) {
    if (score_pair->second.is_decoy()) {
      decoy_score_file << score_pair->first << endl;
      decoy_scores_array[num_decoys++] = score_pair->first;
    } else {
      target_score_file << score_pair->first << endl;
      target_scores_array[num_targets++] = score_pair->first;
    }
  }
  FLOAT_T eta_target = 0.0;
  FLOAT_T beta_target = 0.0;
  FLOAT_T shift_target = 0.0;
  FLOAT_T correlation_target = 0.0;

  FLOAT_T eta_decoy = 0.0;
  FLOAT_T beta_decoy = 0.0;
  FLOAT_T shift_decoy = 0.0;
  FLOAT_T correlation_decoy = 0.0;

  FLOAT_T y;

  // plot fit for targets
  hhc_estimate_weibull_parameters_from_xcorrs(target_scores_array, num_targets, &eta_target, 
	&beta_target, &shift_target, &correlation_target, spectrum, charge);
  for (FLOAT_T x = scores.front().first; x <= scores.back().first; x = x + 0.01) {
      y = (beta_target / eta_target) * pow(((x+shift_target)/eta_target), beta_target - 1) * exp(- pow((x+shift_target)/eta_target, beta_target));
      target_fit_file << x << "\t" << y << endl;
  }

  // plot fit for decoys 
  hhc_estimate_weibull_parameters_from_xcorrs(decoy_scores_array, num_decoys, &eta_decoy, 
	&beta_decoy, &shift_decoy, &correlation_decoy, spectrum, charge);
  for (FLOAT_T x = scores.front().first; x <= scores.back().first; x = x + 0.01) {
      y = (beta_decoy / eta_decoy) * pow(((x+shift_decoy)/eta_decoy), beta_decoy - 1) 
	* exp(- pow((x+shift_decoy)/eta_decoy, beta_decoy));
      decoy_fit_file << x << "\t" << y << endl;
  }

  cout << "target correlation " << correlation_target << " <br>" << endl;
  cout << "decoy correlation " << correlation_decoy << " <br>" << endl;

  FLOAT_T pvalue_target;
  FLOAT_T pvalue_decoy;

  for (vector<pair<FLOAT_T, LinkedPeptide> >::reverse_iterator score_pair = scores.rbegin();
	score_pair != scores.rend(); ++score_pair) {
    pvalue_target = score_logp_weibull(score_pair->first, eta_target, beta_target);
    pvalue_decoy = score_logp_weibull(score_pair->first, eta_decoy, beta_decoy);
    if (score_pair->first+shift_target <= 0) 
      pvalue_target = 1;
    else  
      pvalue_target = exp( - pow( (score_pair->first+shift_target)/eta_target, beta_target));
    if (score_pair->first+shift_decoy <= 0) 
      pvalue_decoy = 1;
    else 
      pvalue_decoy = exp( - pow( (score_pair->first+shift_decoy)/eta_decoy, beta_decoy));
    if (score_pair->second.is_decoy()) {
      target_pvalue_file << pvalue_target << endl;
      decoy_pvalue_file << pvalue_decoy << endl;
    }
  }
} 
/*
string shuffle(string other) {
  string shuffled = string(other);
  int start_idx = 1;
  int end_idx = other.size() - 2;
  int switch_idx = 0;
  char temp_char = 0;
  while(start_idx < end_idx){
    switch_idx = get_random_number_interval(start_idx, end_idx);
    temp_char = shuffled[start_idx];
    shuffled[start_idx] = shuffled[switch_idx];
    shuffled[switch_idx] = temp_char;
    ++start_idx;
  }
  return shuffled;
}
*/

// shuffle a Peptide, preserving N and C terminals
Peptide shuffle(Peptide peptide) {
  string shuffled = string(peptide.sequence());
  Peptide shuffled_peptide = Peptide();
  int start_idx = 1;
  int end_idx = peptide.length() - 2;
  int switch_idx = 0;
  char temp_char = 0;
  while(start_idx < end_idx){
    switch_idx = get_random_number_interval(start_idx, end_idx);
    temp_char = shuffled[start_idx];
    shuffled[start_idx] = shuffled[switch_idx];
    shuffled[switch_idx] = temp_char;
    if (peptide.has_link_at(switch_idx))
      shuffled_peptide.add_link(start_idx, peptide.link_at(switch_idx));
    if (peptide.has_link_at(start_idx))
      shuffled_peptide.add_link(switch_idx, peptide.link_at(start_idx));
    ++start_idx;
  }
  shuffled_peptide.set_sequence(shuffled);
  return shuffled_peptide;
}


void add_decoys(vector<LinkedPeptide>& decoys, LinkedPeptide& lp) {
  vector<Peptide> peptides = lp.peptides();
  LinkedPeptide decoy = LinkedPeptide(lp.charge(), lp.linker_mass());
  Peptide pepA_shuffled = shuffle(peptides[0]);
  decoy.add_peptide(pepA_shuffled);
  if (lp.size() == 2) {
    Peptide pepB_shuffled = shuffle(peptides[1]);
    decoy.add_peptide(pepB_shuffled);
  }
  decoy.set_decoy();
  decoys.push_back(decoy);
}
/*
void add_decoys(vector<LinkedPeptide>& decoys, LinkedPeptide& lp, char* links, int charge, FLOAT_T linker_mass) {
  vector<Peptide> peptides = lp.peptides();
  set<string> shuffled_peptides;
  shuffled_peptides.insert(shuffle(peptides[0].sequence()));
  if (lp.size() == 2)
    shuffled_peptides.insert(shuffle(peptides[1].sequence()));
  //if (lp.size() < 2) cout << "foooo" << endl; 
  add_linked_peptides(decoys, shuffled_peptides, string(links), linker_mass, charge, true);  

}
*/
/*
void add_decoys(vector<LinkedPeptide>& decoys, LinkedPeptide& lp, char* links, int charge, FLOAT_T linker_mass) {
  vector<Peptide> peps = lp.peptides();
  if (lp.size() == 1) return;
  map<char, set<char> > bond_map;
  string bonds_string = string(links);
  for (int i = 0; i < bonds_string.length() - 2; i += 4) {
     bond_map[bonds_string[i]].insert(bonds_string[i+2]);
  }
  string alpha_sequence = shuffle(peps[0].sequence());
  string beta_sequence = shuffle(peps[1].sequence());
  char* sequenceA = (char*) alpha_sequence.c_str();
  char* sequenceB = (char*) beta_sequence.c_str();
  for (int i = 0; i < alpha_sequence.length(); ++i) {
     std::map<char, std::set<char> >::iterator char_it = bond_map.find(alpha_sequence[i]);
     if (char_it != bond_map.end()) {
       for (int j = 0; j < beta_sequence.length(); ++j) {
         if (char_it->second.find(beta_sequence[j]) != char_it->second.end()) {
           //cout << alpha_sequence << " " << i << " " << beta_sequence << " " << j << endl;
	   LinkedPeptide lp = LinkedPeptide(sequenceA,sequenceB,i,j,linker_mass,charge);
	   //cout << "new linked peptide: " << lp << endl;
           decoys.push_back(lp);
          }
         }
        }
      }
   for (int i = 0; i < beta_sequence.length(); ++i) {
     std::map<char, set<char> >::iterator char_it = bond_map.find(beta_sequence[i]);
     if (char_it != bond_map.end()) {
       for (int j = 0; j < alpha_sequence.length(); ++j) {
         if (char_it->second.find(alpha_sequence[j]) != char_it->second.end()) {
	   LinkedPeptide lp = LinkedPeptide(sequenceB,sequenceA,i,j,linker_mass,charge);
           decoys.push_back(lp);
         }
       }	
     }
    }
}*/

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
        //cout << "peptide " << get_peptide_sequence(peptide) << endl;
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
			     int charge,
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
  string bonds_string = string(links);
  set<string> peptides;
  get_linkable_peptides(peptides, protein_iterator, peptide_constraint);
  add_linked_peptides(all_ions, peptides, bonds_string, linker_mass, charge, false); 
}


void add_linked_peptides(vector<LinkedPeptide>& all_ions, set<string>& peptides, string links, FLOAT_T linker_mass, int charge, bool is_decoy) {
  //if (peptides.size() < 2) return;
  bool single_decoy = false;
  BondMap bonds; 
  vector<LinkedPeptide> ions;
  for (int i = 0; i < links.length() - 2; i += 4) {
     bonds[links[i+2]].insert(links[i]);
     bonds[links[i]].insert(links[i+2]);
  }
  //if (is_decoy) cout << "adding decoy" << endl; else cout << "adding real" << endl;  
  if (peptides.size() == 1) {
    single_decoy = true;
    carp(CARP_DETAILED_DEBUG, "single decoy");
  }
  // iterate over both sequences, adding linked peptides with correct links
  for (set<string>::iterator pepA = peptides.begin(); pepA != peptides.end(); ++pepA) {
    char* sequenceA = (char*) pepA->c_str();
    // add unlinked precursor
    LinkedPeptide lp = LinkedPeptide(charge, linker_mass);
    Peptide p = Peptide(sequenceA);
    lp.add_peptide(p);
    if (!is_decoy || single_decoy) ions.push_back(lp);
    for (int i = 0; i < pepA->length(); ++i) {
      BondMap::iterator bond = bonds.find(pepA->at(i));
      // if a link aa and doesn't end in K
      if (bond != bonds.end() && i != pepA->length()-1) {
	if (i == pepA->length()-1 && pepA->at(pepA->length()-1) == 'K') continue;
        // add dead end
	if (!is_decoy || single_decoy) ions.push_back(LinkedPeptide(sequenceA, NULL, i, -1, linker_mass, charge));
        // add self loop
	for (int j = i+1; j < pepA->length(); ++j) {
          if (bond->second.find(pepA->at(j)) != bond->second.end()) { 
	    //skip if linked to a K at the end
	    if (j == pepA->length()-1 && pepA->at(pepA->length()-1) == 'K') continue;
	    if (!is_decoy || single_decoy) ions.push_back(LinkedPeptide(sequenceA, NULL, i, j, linker_mass, charge));
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
	      if (!is_decoy || strcmp(sequenceA, sequenceB) != 0)
	      ions.push_back(LinkedPeptide(sequenceA, sequenceB, i, j, linker_mass, charge));
	    }
          }
	} // get next pepB 
      }
    }
  } // get next pepA 
  for (vector<LinkedPeptide>::iterator ion = ions.begin(); ion != ions.end(); ++ion) {
    if (is_decoy) { 
      ion->set_decoy();
      //cout << "decoy\t" << *ion << endl;
    } else {
      //cout << "real\t" << *ion << endl;
    }
   }
   all_ions.insert(all_ions.end(), ions.begin(), ions.end());
}

BOOLEAN_T hhc_estimate_weibull_parameters_from_xcorrs(
  FLOAT_T* scores,
  int num_scores,
  FLOAT_T* eta,
  FLOAT_T* beta,
  FLOAT_T* shift,
  FLOAT_T* correlation,
  SPECTRUM_T* spectrum,
  int charge
  ){
/*
  if( match_collection == NULL || spectrum == NULL ){
    carp(CARP_ERROR, "Cannot estimate parameters from null inputs.");
    return FALSE;
  }
*/
  // check that we have the minimum number of matches
  if( num_scores < MIN_WEIBULL_MATCHES ){
    carp(CARP_DETAILED_INFO, "Too few psms (%i) to estimate "
         "p-value parameters for spectrum %i, charge %i",
         num_scores, get_spectrum_first_scan(spectrum), charge);
    // set eta, beta, and shift to something???
    return FALSE;
  }

  // randomly sample n from the list by shuffling and taking first n 
  shuffle_floats(scores, num_scores);
  int num_samples = num_scores;
  //int num_samples = PARAM_ESTIMATION_SAMPLE_COUNT;
  //if(num_samples > num_scores){ num_samples = num_scores; }

  //int num_samples = num_scores;
  // reverse sort the first num_samples of them
  qsort(scores, num_samples, sizeof(float), compare_floats_descending);

  // use only a fraction of the samples, the high-scoring tail
  // this parameter is hidden from the user
  double fraction_to_fit = get_double_parameter("fraction-top-scores-to-fit");
  assert( fraction_to_fit >= 0 && fraction_to_fit <= 1 );
  int num_tail_samples = (int)(num_samples * fraction_to_fit);
  carp(CARP_DEBUG, "Estimating Weibull params with %d psms (%.2f of %i)",
       num_tail_samples, fraction_to_fit, num_samples);

  // do the estimation
  fit_three_parameter_weibull(scores, num_tail_samples, num_samples,
      MIN_XCORR_SHIFT, MAX_XCORR_SHIFT, XCORR_SHIFT,
      eta, beta, shift, correlation);

  carp(CARP_DETAILED_DEBUG,
      "Correlation: %.6f\nEta: %.6f\nBeta: %.6f\nShift: %.6f\n",
      correlation, *eta, *beta, *shift);

  return TRUE;
}
