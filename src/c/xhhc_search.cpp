#include "hhc_ion_series.h"
#include "objects.h"
#include <fstream>
#include <math.h>
//#include <iostream>

// get rid of these
#define PARAM_ESTIMATION_SAMPLE_COUNT 500
#define MIN_WEIBULL_MATCHES 40 
#define MIN_XCORR_SHIFT -5.0
#define MAX_XCORR_SHIFT  5.0
#define XCORR_SHIFT 0.05

// mine
#define BONFERRONI_CUT_OFF_P 0.0001
#define BONFERRONI_CUT_OFF_NP 0.01
#define MIN_WEIBULL_SAMPLES 1000
#define MIN_PRECURSORS 3
using namespace std;

typedef map<char, set<char> > BondMap;

double bonf_correct_2(double p_value);

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

void add_linked_peptides(vector<LinkedPeptide>& all_ions, set<string>& peptides, string links, FLOAT_T linker_mass, int charge);

void get_ions_from_mz_range(vector<LinkedPeptide>& filtered_ions,
	vector<LinkedPeptide>& all_ions,
	FLOAT_T precursor_mass,
	int charge,
	int mass_window,
	int decoy_iterations);

void add_decoys(vector<LinkedPeptide>& decoys, LinkedPeptide& lp);

void plot_weibull(vector<pair<FLOAT_T, LinkedPeptide> >& scores, SPECTRUM_T* spectrum, int charge); 

void find_all_precursor_ions(vector<LinkedPeptide>& all_ions, 
	char* links, 
	FLOAT_T linker_mass,
	char* missed_link_cleavage, 
	char* database_file,
	int charge);

int main(int argc, char** argv) {
  char* missed_link_cleavage = "K";
  int num_missed_cleavages = 0;
  char* ms2_file = NULL;
  //char* min_mass_string = NULL;
  //char* max_mass_string = NULL;
  int mass_window = 2;
  char* database = NULL;
  char* links = NULL;
  char* linker_mass_string = NULL;
  int decoy_iterations = 5;
  int scan_num = 0;
  int charge = 1;
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

  parse_arguments_set_req(
    "linker mass", 
    "combined mass of linker and linker modifications", 
    (void *) &linker_mass_string, 
    STRING_ARG);

  parse_arguments_set_req(
    "ms2-filename", 
    "A file containing multiple MS-MS spectra in .ms2 format.",
    (void *) &ms2_file,
    STRING_ARG);
 
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
    "mass-window",
    "window around precursor m/z to select target peptides", 
    (void *) &mass_window,
    INT_ARG); 

  initialize_parameters();
  parse_arguments(argc, argv, 0);
  // something wrong with DOUBLE_ARG
  FLOAT_T linker_mass = atof(linker_mass_string);

  vector<LinkedPeptide> all_ions;
  
  find_all_precursor_ions(all_ions, links, linker_mass, missed_link_cleavage, database,1);

  // sort filtered ions and decoy ions by mass
  sort(all_ions.begin(), all_ions.end());

  //SPECTRUM_T* spectrum = allocate_spectrum();
  SPECTRUM_COLLECTION_T* spectra = new_spectrum_collection(ms2_file);
  parse_spectrum_collection(spectra);
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* spectrum_iterator = 
	new_filtered_spectrum_charge_iterator(spectra);
 
  FLOAT_T score;
 // best pvalues
  double pvalue_decoy;
  double pvalue_target;

  double pvalue_target_single;
  double pvalue_decoy_single;
  double pvalue_target_linked;
  double pvalue_decoy_linked;
 
  ofstream search_target_file ("crux-output3/search.target.txt");
  ofstream search_decoy_file ("crux-output3/search.decoy.txt");
  // main loop over spectra in ms2 file
  while (filtered_spectrum_charge_iterator_has_next(spectrum_iterator)) {
    SPECTRUM_T* spectrum = filtered_spectrum_charge_iterator_next(spectrum_iterator, &charge);
    SCORER_T* scorer = new_scorer(XCORR);
    scan_num = get_spectrum_first_scan(spectrum);

    cout << "scan " << scan_num << endl;
    
    //vector<pair<FLOAT_T, LinkedPeptide> > linked_scores;
    //vector<pair<FLOAT_T, LinkedPeptide> > single_scores;
    vector<pair<FLOAT_T, LinkedPeptide> > scores;

    vector<LinkedPeptide> filtered_ions;

    FLOAT_T precursor_mass = get_spectrum_neutral_mass(spectrum, charge); 

    get_ions_from_mz_range(
	filtered_ions, // stored in this vector
	all_ions,
	precursor_mass,
	charge,
	mass_window,
	decoy_iterations);
    if (filtered_ions.size() < MIN_PRECURSORS) {
      cout << "not enough precursors found in range, skipping" << endl;
      continue;
    }
    
    cout << "num ions " << filtered_ions.size() << endl;
    LinkedIonSeries ion_series;
    // keep track of number of decoys for estimating pvalue  
    int num_decoys = 0;
    int num_single_decoys = 0;
    int num_linked_decoys = 0;
    /// for every ion in the mass window
    for (vector<LinkedPeptide>::iterator ion = filtered_ions.begin();
      ion != filtered_ions.end(); ++ion) {
      ion_series = LinkedIonSeries(links, charge, linker_mass);
      ion_series.add_linked_ions(*ion);
      score = hhc_score_spectrum_v_ion_series(scorer, spectrum, ion_series);
      //cout << "score for " << *ion << " is " << score << endl;

      if (ion->is_decoy()) {
        ++num_decoys;
        if (ion->is_single())
          ++num_single_decoys;
        else
          ++num_linked_decoys;
      }

      scores.push_back(make_pair(score, *ion));
      //cout << "size of scores is now " << scores.size() << endl;
    }
    
    //sort(linked_scores.begin(), linked_scores.end());
    //sort(single_scores.begin(), single_scores.end());
    sort(scores.begin(), scores.end());

    cout << "num decoys " << num_decoys << endl;
    cout << "num_single_decoys " << num_single_decoys << endl;
    cout << "num_linked_decoys " << num_linked_decoys << endl;
    // add enough decoys to estimate pvalues
    
    vector<LinkedPeptide> extra_single_decoys;
    vector<LinkedPeptide> extra_linked_decoys;
    //vector<LinkedPeptide> extra_decoys;

    if (num_single_decoys == 0) {
	num_single_decoys = MIN_WEIBULL_SAMPLES;
	cout << "no single decoys" << endl;
    }
    if (num_linked_decoys == 0)  {
	num_linked_decoys = MIN_WEIBULL_SAMPLES;
	cout << "no linked decoys" << endl;
    }
  
    while (num_single_decoys < MIN_WEIBULL_SAMPLES || num_linked_decoys < MIN_WEIBULL_SAMPLES) {
      //cout << "num_decoys now " << num_decoys << endl;
      for (vector<LinkedPeptide>::iterator peptide = filtered_ions.begin();
	peptide != filtered_ions.end(); ++peptide) {
	  if (!peptide->is_decoy()) { 
            if (num_single_decoys < MIN_WEIBULL_SAMPLES) {
              add_decoys(extra_single_decoys, *peptide);
              ++num_single_decoys;
            }
            if (num_linked_decoys < MIN_WEIBULL_SAMPLES) {
              add_decoys(extra_linked_decoys, *peptide);
              ++num_linked_decoys;
            }
	    //add_decoys(extra_decoys, *peptide);
	    //++num_decoys;
          }
	} 
    }
    cout << "num single decoy scores " << num_single_decoys << endl;
    cout << "num linked decoy scores " << num_linked_decoys << endl;
    //cout << "extra " << extra_decoys.size() << endl;
    //FLOAT_T decoy_scores_array[num_decoys];
    FLOAT_T single_decoy_scores_array[num_single_decoys];
    FLOAT_T linked_decoy_scores_array[num_linked_decoys];

    int single_decoy_index = 0;
    int linked_decoy_index = 0;

    for (int i = 0; i < scores.size(); ++i) { 
      LinkedPeptide& lp = scores[i].second;
      if (lp.is_decoy()) {
	if (lp.is_single()) {
          single_decoy_scores_array[single_decoy_index++] = scores[i].first;
        } else {
          linked_decoy_scores_array[linked_decoy_index++] = scores[i].first;
        }
        //decoy_scores_array[decoy_index++] = scores[i].first;
      }
    }
    cout << "single_decoy_index " << single_decoy_index << endl; 
    cout << "linked_decoy_index " << linked_decoy_index << endl; 
    for (vector<LinkedPeptide>::iterator extra_decoy = extra_single_decoys.begin(); extra_decoy != extra_single_decoys.end(); ++extra_decoy) {
      ion_series = LinkedIonSeries(links,charge, linker_mass);
      ion_series.add_linked_ions(*extra_decoy);
      score = hhc_score_spectrum_v_ion_series(scorer, spectrum, ion_series);
      single_decoy_scores_array[single_decoy_index++] = score;
    }

    for (vector<LinkedPeptide>::iterator extra_decoy = extra_linked_decoys.begin(); extra_decoy != extra_linked_decoys.end(); ++extra_decoy) {
      ion_series = LinkedIonSeries(links,charge, linker_mass);
      ion_series.add_linked_ions(*extra_decoy);
      score = hhc_score_spectrum_v_ion_series(scorer, spectrum, ion_series);
      linked_decoy_scores_array[linked_decoy_index++] = score;
    }

    // find top single and linked target and decoy psm
    vector<pair<FLOAT_T, LinkedPeptide> >::reverse_iterator best_single_target = scores.rbegin();
    vector<pair<FLOAT_T, LinkedPeptide> >::reverse_iterator best_single_decoy = scores.rbegin();
    vector<pair<FLOAT_T, LinkedPeptide> >::reverse_iterator best_linked_target = scores.rbegin();
    vector<pair<FLOAT_T, LinkedPeptide> >::reverse_iterator best_linked_decoy = scores.rbegin();

    vector<pair<FLOAT_T, LinkedPeptide> >::reverse_iterator best_target;
    vector<pair<FLOAT_T, LinkedPeptide> >::reverse_iterator best_decoy;
    if (single_decoy_index != 0) {
      while (best_single_target->second.is_decoy() || !best_single_target->second.is_single()) 
        ++best_single_target;
      while (!best_single_decoy->second.is_decoy() || !best_single_target->second.is_single()) 
        ++best_single_decoy;
    }
    if (linked_decoy_index != 0) {
      while (best_linked_target->second.is_decoy() || best_linked_target->second.is_single()) 
        ++best_linked_target;
      while (!best_linked_decoy->second.is_decoy() || best_linked_target->second.is_single()) 
        ++best_linked_decoy;
    }
   // weibull parameters for single and cross-linked 
    FLOAT_T eta_single = 0.0;
    FLOAT_T beta_single = 0.0;
    FLOAT_T shift_single  = 0.0;
    FLOAT_T correlation_single  = 0.0;

    FLOAT_T eta_linked = 0.0;
    FLOAT_T beta_linked  = 0.0;
    FLOAT_T shift_linked  = 0.0;
    FLOAT_T correlation_linked  = 0.0;
    // fit weibull to decoys
    hhc_estimate_weibull_parameters_from_xcorrs(single_decoy_scores_array, num_single_decoys, &eta_single, 
	&beta_single, &shift_single, &correlation_single, spectrum, charge);

    hhc_estimate_weibull_parameters_from_xcorrs(linked_decoy_scores_array, num_linked_decoys, &eta_linked, 
	&beta_linked, &shift_linked, &correlation_linked, spectrum, charge);

    // get p-values
    if (single_decoy_index == 0) {
      pvalue_target_single = -1; 
      pvalue_decoy_single = -1;
    } else {
      pvalue_target_single = bonf_correct_2(score_logp_bonf_weibull(best_single_target->first,
	eta_single, beta_single, shift_single, filtered_ions.size()));
      pvalue_decoy_single  = bonf_correct_2(score_logp_bonf_weibull(best_single_decoy->first, 
	eta_single, beta_single, shift_single, filtered_ions.size()));
    }

    if (linked_decoy_index == 0) {
      pvalue_target_linked = -1;
      pvalue_decoy_linked = -1;
    } else {
      pvalue_target_linked = bonf_correct_2(score_logp_bonf_weibull(best_linked_target->first, 
	eta_linked, beta_linked, shift_linked, filtered_ions.size()));
      pvalue_decoy_linked  = bonf_correct_2(score_logp_bonf_weibull(best_linked_decoy->first, 
	eta_linked, beta_linked, shift_linked, filtered_ions.size()));
    }

    // find best overall target and decoy
    if (pvalue_target_single > pvalue_target_linked) {
      cout << "pvalue target single > pvalue target linked " << endl;
      best_target = best_single_target;
      pvalue_target = pvalue_target_single;
    } else {
      cout << "pvalue target single < pvalue target linked " << endl;
      best_target = best_linked_target;
    cout << "best_linked_target " << best_linked_target->first << "\t" << best_linked_target->second << endl;
      cout << "\t\t" << best_target->second << endl; 
      pvalue_target = pvalue_target_linked;
    } 
    if (pvalue_decoy_single > pvalue_decoy_linked) {
      best_decoy = best_single_decoy;
      pvalue_decoy = pvalue_decoy_single;
    } else {
      best_decoy = best_linked_decoy;
      pvalue_decoy = pvalue_decoy_linked;
    }
    cout << "pvalue_target " << pvalue_target << endl;
    cout << "pvalue_decoy " << pvalue_decoy << endl;
  /*  if (best_target->first+shift <= 0) 
      pvalue_target = 1;
    else  
      pvalue_target = exp( - pow( (best_target->first+shift)/eta, beta));
    if (best_decoy->first+shift <= 0)
      pvalue_decoy = 1;
    else
      pvalue_decoy = exp( - pow( (best_decoy->first+shift)/eta, beta));
    */
    cout << best_target->second << endl;
    cout << "T  xcorr: " << best_target->first << "\t-log(pvalue): " << 
	pvalue_target << "\t" << best_target->second << endl;
    cout << "D  xcorr: " << best_decoy->first << "\t-log(pvalue): " << 
	pvalue_decoy << "\t" << best_decoy->second << endl;
    
    search_target_file << 
	scan_num << "\t" <<
	charge << "\t" <<
	best_target->second.get_mz() << "\t" <<
	best_target->second << "\t" <<
	best_target->first << "\t" <<
	pvalue_target << endl;

    search_decoy_file << 
	scan_num << "\t" <<
	charge << "\t" <<
	best_decoy->second.get_mz() << "\t" <<
	best_decoy->second << "\t" <<
	best_decoy->first << "\t" <<
	pvalue_decoy << endl;
    //free_scorer(scorer);
    //free_spectrum(spectrum);
    
  }
  search_target_file.close();
  search_decoy_file.close();
  //free_scorer(scorer);
  //free_spectrum_collection(spectra);
  //free_spectrum(spectrum);
  return 0;
}

double bonf_correct_2(double pvalue) {
  double real_pvalue = exp(-pvalue);
  if(real_pvalue > BONFERRONI_CUT_OFF_P 
     || real_pvalue*2 > BONFERRONI_CUT_OFF_NP){

    double corrected_pvalue = -log(1-pow((1-real_pvalue), 2));
    carp(CARP_DETAILED_DEBUG, "Stat: pvalue after = %.6f", corrected_pvalue);
    return corrected_pvalue;
  }
  // else, use the approximation
  else{
    double corrected_pvalue = -log(real_pvalue*2);
    carp(CARP_DETAILED_DEBUG, "Stat: pvalue after = %.6f", corrected_pvalue);
    return corrected_pvalue;
  }
}
void get_ions_from_mz_range(vector<LinkedPeptide>& filtered_ions,
	vector<LinkedPeptide>& all_ions,
	FLOAT_T precursor_mass,
	int charge,
	int mass_window,
	int decoy_iterations) {
  FLOAT_T min_mass = precursor_mass - mass_window;
  FLOAT_T max_mass = precursor_mass + mass_window;
  cout << "min mass " << min_mass << " max mass " << max_mass << endl;
  FLOAT_T ion_mass;
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin();
	ion != all_ions.end(); ++ion) {
    ion->set_charge(charge);
    ion->calculate_mass();
    ion_mass = ion->mass();
    if (ion_mass >= min_mass && ion_mass <= max_mass) {
      filtered_ions.push_back(*ion);
      //cout << *ion << "\t" << ion->mass() << endl;
      for (int i = decoy_iterations; i > 0; --i)
        add_decoys(filtered_ions, *ion);
    }
  }
}

Peptide shuffle(Peptide peptide) {
  Peptide p = Peptide();
  string shuffled = string(peptide.sequence());
  //cout << shuffled << " before shuffle " << endl;
  Peptide shuffled_peptide = Peptide();
  for (int i = 0; i < shuffled.length(); ++i) {
    if (peptide.has_link_at(i)) shuffled_peptide.add_link(i, p);
  }

  int start_idx = 1;
  int end_idx = peptide.length() - 2;
  int switch_idx = 0;
  char temp_char = 0;
  while(start_idx <= end_idx){
    //cout << shuffled << endl;
    switch_idx = get_random_number_interval(start_idx, end_idx);
    //cout << "start index " << start_idx << " switch index " << switch_idx << endl;
    temp_char = shuffled[start_idx];
    shuffled[start_idx] = shuffled[switch_idx];
    shuffled[switch_idx] = temp_char;
    if (shuffled_peptide.has_link_at(switch_idx)) {
      //if not a self loop
      if (!shuffled_peptide.has_link_at(start_idx)) {
        //cout << "there is link at " << switch_idx << " to a " << shuffled[switch_idx] << endl;
        shuffled_peptide.remove_link(switch_idx);
        shuffled_peptide.add_link(start_idx, p);
      }
    } else if (shuffled_peptide.has_link_at(start_idx)) {
      //cout << "there is link at " << start_idx << " to a " << shuffled[start_idx] << endl;
      shuffled_peptide.remove_link(start_idx);
      shuffled_peptide.add_link(switch_idx, p);
    }
    ++start_idx;
  }
  //cout << shuffled << " after shuffle " << endl;
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
  //cout << "decoy " << decoy << " from " << endl << "peptd " << lp << endl;
  //cout << "decoy\ttarget" << endl;
  //cout << decoy.get_mz() << "\t" << lp.get_mz() << endl;
  //cout << decoy.mass() << "\t" << lp.mass() << endl;
  decoy.set_decoy();
  decoy.calculate_mass();
  decoys.push_back(decoy);
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
			     char* missed_link_cleavage,
		             char* database_file,
			     int charge)
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
  add_linked_peptides(all_ions, peptides, bonds_string, linker_mass, charge); 
}


void add_linked_peptides(vector<LinkedPeptide>& all_ions, set<string>& peptides, string links, FLOAT_T linker_mass, int charge) {
  //if (peptides.size() < 2) return;
  BondMap bonds; 
  vector<LinkedPeptide> ions;
  for (int i = 0; i < links.length() - 2; i += 4) {
     bonds[links[i+2]].insert(links[i]);
     bonds[links[i]].insert(links[i+2]);
  }

  // iterate over both sequences, adding linked peptides with correct links
  for (set<string>::iterator pepA = peptides.begin(); pepA != peptides.end(); ++pepA) {
    char* sequenceA = (char*) pepA->c_str();
    // add unlinked precursor
    LinkedPeptide lp = LinkedPeptide(charge, linker_mass);
    Peptide p = Peptide(sequenceA);
    lp.add_peptide(p);
    ions.push_back(lp);
    for (int i = 0; i < pepA->length(); ++i) {
      BondMap::iterator bond = bonds.find(pepA->at(i));
      // if a link aa and doesn't end in K
      if (bond != bonds.end() && i != pepA->length()-1) {
	if (i == pepA->length()-1 && pepA->at(pepA->length()-1) == 'K') continue;
        // add dead end
	ions.push_back(LinkedPeptide(sequenceA, NULL, i, -1, linker_mass, charge));
        // add self loop
	for (int j = i+1; j < pepA->length(); ++j) {
          if (bond->second.find(pepA->at(j)) != bond->second.end()) { 
	    //skip if linked to a K at the end
	    if (j == pepA->length()-1 && pepA->at(pepA->length()-1) == 'K') continue;
	    ions.push_back(LinkedPeptide(sequenceA, NULL, i, j, linker_mass, charge));
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
	      if (strcmp(sequenceA, sequenceB) != 0)
	      ions.push_back(LinkedPeptide(sequenceA, sequenceB, i, j, linker_mass, charge));
	    }
          }
	} // get next pepB 
      }
    }
  } // get next pepA 
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
