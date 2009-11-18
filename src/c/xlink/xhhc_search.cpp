#include "xhhc_scorer.h"
#include "xhhc_ion_series.h"
//#include "xhhc_search.h"

//CRUX INCLUDES
extern "C" {
#include "objects.h"
}


#include <fstream>
#include <math.h>
//#include <iostream>
/*
// get rid of these
#define PARAM_ESTIMATION_SAMPLE_COUNT 500
#define MIN_WEIBULL_MATCHES 40 
#define MIN_XCORR_SHIFT -5.0
#define MAX_XCORR_SHIFT  5.0
#define XCORR_SHIFT 0.05

// mine
#define BONFERRONI_CUT_OFF_P 0.0001
#define BONFERRONI_CUT_OFF_NP 0.01
#define MIN_WEIBULL_SAMPLES 750 
#define MIN_PRECURSORS 3
*/
using namespace std;

typedef map<char, set<char> > BondMap;

double bonf_correct(double nlp_value, int nt);




void get_ions_from_mz_range(vector<LinkedPeptide>& filtered_ions,
	vector<LinkedPeptide>& all_ions,
	FLOAT_T precursor_mass,
	int charge,
	FLOAT_T mass_window,
	int decoy_iterations);

void plot_weibull(vector<pair<FLOAT_T, LinkedPeptide> >& scores, SPECTRUM_T* spectrum, int charge); 

int main(int argc, char** argv) {
  char* missed_link_cleavage = "K";
  int num_missed_cleavages = 0;
  char* ms2_file = NULL;
  //char* min_mass_string = NULL;
  //char* max_mass_string = NULL;

  char* str_mass_window = NULL;
  char* str_mass_window_decoy = NULL;

  FLOAT_T mass_window = 5;
  FLOAT_T mass_window_decoy = 5;

  char* database = NULL;
  char* links = NULL;
  char* linker_mass_string = NULL;
  int decoy_iterations = 5;
  int min_weibull_points = 400;
  int scan_num = 0;
  int charge = 1;

  int top_match = 1;

  //bool open_modification = false;
  BOOLEAN_T compute_pvalues = FALSE;
  //int open_modification_int = 0;
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
    "compute-p-values",
    "",
    (void *) &compute_pvalues,
    BOOLEAN_ARG);
 
  parse_arguments_set_opt(
    "decoy-iterations",
    "",
    (void *) &decoy_iterations,
    INT_ARG);

  parse_arguments_set_opt(
    "min-weibull-points",
    "",
    (void *)& min_weibull_points,
    INT_ARG);



  parse_arguments_set_opt(
    "missed-link-cleavage",
    "",
    (void *) &missed_link_cleavage, 
    STRING_ARG);

  parse_arguments_set_opt(
    "top-match",
    "number of psms to print per scan and charge",
    (void*)&top_match,
    INT_ARG);

  parse_arguments_set_opt(
    "num-missed-cleavages", 
    "maximum number of missed cleavages (not including one at link site)", 
    (void *) &num_missed_cleavages, 
    INT_ARG);


  parse_arguments_set_opt(
    "mass-window",
    "window around precursor m/z to select target peptides", 
    (void *) &str_mass_window,
    STRING_ARG); 

  parse_arguments_set_opt(
    "mass-window-decoy",
    "window around precursor m/z to select decoy peptides",
    (void *) &str_mass_window_decoy,
    STRING_ARG);

  initialize_parameters();
  if (!parse_arguments(argc, argv, 0)) {
   char* error_message;
   char* usage = parse_arguments_get_usage("xhhc-search");
   int result = parse_arguments_get_error(&error_message);
   fprintf(stderr, "Error in command line. Error # %d\n", result);
   fprintf(stderr, "%s\n", error_message);
   fprintf(stderr, "%s", usage);
   free(usage);
   exit(1);
 }
  // something wrong with DOUBLE_ARG
  FLOAT_T linker_mass = atof(linker_mass_string);

  mass_window = atof(str_mass_window);
  mass_window_decoy = atof(str_mass_window_decoy);
  
  LinkedPeptide::linker_mass = linker_mass;
  vector<LinkedPeptide> all_ions;
  
  find_all_precursor_ions(all_ions, links, missed_link_cleavage, database,1);

  // sort filtered ions and decoy ions by mass
  sort(all_ions.begin(), all_ions.end());

  SPECTRUM_T* spectrum = allocate_spectrum();
  SPECTRUM_COLLECTION_T* spectra = new_spectrum_collection(ms2_file);
  parse_spectrum_collection(spectra);
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* spectrum_iterator = 
	new_filtered_spectrum_charge_iterator(spectra);
 
  FLOAT_T score;
 // best pvalues

 
  ofstream search_target_file ("crux-output/search.target.txt");
  //print header
  search_target_file << "scan\t";
  search_target_file << "charge\t";
  search_target_file << "spectrum precursor m/z\t";
  search_target_file << "spectrum neutral mass\t";
  search_target_file << "peptide mass\t";
  search_target_file << "xcorr score\t";
  search_target_file << "xcorr rank\t";
  search_target_file << "p-value\t";
  search_target_file << "matches/spectrum\t";
  search_target_file << "sequence"<<endl;
  


  ofstream search_decoy_file ("crux-output/search.decoy.txt");
  //print header
  search_decoy_file << "scan\t";
  search_decoy_file << "charge\t";
  search_decoy_file << "spectrum precursor m/z\t";
  search_decoy_file << "spectrum neutral mass\t";
  search_decoy_file << "peptide mass\t";
  search_decoy_file << "xcorr score\t";
  search_decoy_file << "xcorr rank\t";
  search_decoy_file << "p-value\t";
  search_decoy_file << "matches/spectrum\t";
  search_decoy_file << "sequence"<<endl;

  Scorer hhc_scorer;
  // main loop over spectra in ms2 file
 
  // for every observed spectrum 
  while (filtered_spectrum_charge_iterator_has_next(spectrum_iterator)) {
    spectrum = filtered_spectrum_charge_iterator_next(spectrum_iterator, &charge);
    //SCORER_T* scorer = new_scorer(XCORR);
    scan_num = get_spectrum_first_scan(spectrum);

    cout << "scan " << scan_num << endl;
    
    //vector<pair<FLOAT_T, LinkedPeptide> > linked_scores;
    //vector<pair<FLOAT_T, LinkedPeptide> > single_scores;
    vector<pair<FLOAT_T, LinkedPeptide> > scores;

    vector<LinkedPeptide> target_xpeptides;
    vector<LinkedPeptide> target_decoy_xpeptides;
    vector<LinkedPeptide> decoy_train_xpeptides;
    vector<LinkedPeptide> decoy_xpeptides;

    FLOAT_T precursor_mz = get_spectrum_precursor_mz(spectrum);
    FLOAT_T precursor_mass = get_spectrum_neutral_mass(spectrum, charge); 




    cout << "finding target xpeptides in mass window..."<<mass_window<<endl;
    get_ions_from_mz_range(
	target_xpeptides, // stored in this vector
	all_ions,
	precursor_mass,
	charge,
	mass_window,
	0);

    if (target_xpeptides.size() < 1) {
      cout << "not enough precursors found in range, skipping" << endl;
      continue;
    }
    

    cout <<"finding training xpeptides in decoy mass window.."<<mass_window_decoy<<endl;
    get_ions_from_mz_range(
	target_decoy_xpeptides,
	all_ions,
	precursor_mass,
	charge,
	mass_window_decoy,
	0);

    cout <<"Creating decoys for target window"<<endl;
    //create the decoys from the target found in the target_mass_window.
    for (vector<LinkedPeptide>::iterator ion = target_xpeptides.begin();
	 ion != target_xpeptides.end(); ++ion) {
        add_decoys(decoy_xpeptides, *ion);
    }
    
    
    cout <<"Creating decoys for decoy mass window"<<endl;
    //create the decoys from the target found in the decoy_mass_window.
    while (decoy_train_xpeptides.size() < min_weibull_points) {
      cout<<"np:"<<decoy_train_xpeptides.size()<<endl;
      for (vector<LinkedPeptide>::iterator ion = target_decoy_xpeptides.begin();
	   ion != target_decoy_xpeptides.end(); ++ion) {
	add_decoys(decoy_train_xpeptides, *ion);
      }
    }    
    
    cout <<"num targets:"<<target_xpeptides.size()<<endl;
    cout <<"num decoys:"<<decoy_xpeptides.size()<<endl;
    cout <<"num training decoys:"<<decoy_train_xpeptides.size()<<endl;



    // for every ion in the mass window
    

    cout <<"Scoring targets"<<endl;
    for (int i=0;i<target_xpeptides.size();i++) {
      LinkedIonSeries ion_series = LinkedIonSeries(links, charge);
      ion_series.add_linked_ions(target_xpeptides[i]);
      score = hhc_scorer.score_spectrum_vs_series(spectrum, ion_series);
      scores.push_back(make_pair(score, target_xpeptides[i]));
    }
    cout <<"Scoring decoys."<<endl;
    for (int i=0;i<decoy_xpeptides.size();i++) {
      LinkedIonSeries ion_series = LinkedIonSeries(links, charge);
      ion_series.add_linked_ions(decoy_xpeptides[i]);
      score = hhc_scorer.score_spectrum_vs_series(spectrum, ion_series);
      scores.push_back(make_pair(score, decoy_xpeptides[i]));
    }

    // add enough decoys to estimate pvalues
    cout << "done" << endl;

    // create arrays to pass to crux's weibull methods

    //FLOAT_T decoy_scores_array[num_decoys];
    FLOAT_T* linked_decoy_scores_array = new FLOAT_T[decoy_train_xpeptides.size()+target_xpeptides.size()];
    
    //use the decoy scores to build the estimator.
    cout << "making arrays...";

    //cout << "single_decoy_index " << single_decoy_index << endl; 
    //cout << "linked_decoy_index " << linked_decoy_index << endl; 
    cout << "scoring training decoys..."<<endl;
    // score all training decoys
    for (int i=0;i<decoy_train_xpeptides.size();i++) {
      LinkedIonSeries ion_series = LinkedIonSeries(links, charge);
      //ion_series.clear();
      ion_series.add_linked_ions(decoy_train_xpeptides[i]);
      score = hhc_scorer.score_spectrum_vs_series(spectrum, ion_series);
      linked_decoy_scores_array[i] = score;
    }
    
    for (int i=0;i<scores.size();i++) {
      if (!scores[i].second.is_decoy())
	linked_decoy_scores_array[i+decoy_train_xpeptides.size()] = scores[i].first;
    }
    

    cout << "done" << endl;
    // find top single and linked target and decoy psm

    cout <<"Sorting scores"<<endl;

    sort(scores.begin(), scores.end());


    cout <<"done"<<endl;
   // weibull parameters for single and cross-linked 
    FLOAT_T eta_linked = 0.0;
    FLOAT_T beta_linked  = 0.0;
    FLOAT_T shift_linked  = 0.0;
    FLOAT_T correlation_linked  = 0.0;

    // fit weibull to decoys

    cout <<"Fittng weibull"<<endl;

    hhc_estimate_weibull_parameters_from_xcorrs(linked_decoy_scores_array, 
						decoy_train_xpeptides.size(), 
						&eta_linked, &beta_linked, 
						&shift_linked, &correlation_linked, 
						spectrum, charge);

    cout <<"done"<<endl;

    int ndecoys = 0;
    int ntargets = 0;
    int score_index = 0;

    while (score_index < scores.size() && (ndecoys < top_match || ntargets < top_match)) {
      if (scores[score_index].second.is_decoy() && ndecoys < top_match) {
	double pvalue = score_logp_bonf_weibull(scores[score_index].first, eta_linked, beta_linked, shift_linked, 1);
	double pvalue_bonf = pvalue;//bonf_correct(pvalue, decoy_xpeptides.size());
	
	if (pvalue != pvalue) {
	  pvalue = 0;
	  pvalue_bonf = 0;
	} else if (pvalue_bonf != pvalue_bonf) {
	  pvalue_bonf = 0;
	}
	
	search_decoy_file << scan_num << "\t"; 
	search_decoy_file << charge << "\t"; 
	search_decoy_file << precursor_mz << "\t";
	search_decoy_file << precursor_mass << "\t";
	search_decoy_file << scores[score_index].second.mass() << "\t";
	search_decoy_file << scores[score_index].first <<"\t";
	search_decoy_file << (ndecoys+1) << "\t";
	search_decoy_file << pvalue << "\t";
	search_decoy_file << decoy_xpeptides.size() << "\t";
	search_decoy_file << scores[score_index].second<<endl;

	ndecoys++;
      } else if (ntargets < top_match) {
	ntargets++;

	double pvalue = score_logp_bonf_weibull(scores[score_index].first, eta_linked, beta_linked, shift_linked, 1);
	double pvalue_bonf = pvalue;//bonf_correct(pvalue, target_xpeptides.size());
	
	if (pvalue != pvalue) {
	  pvalue = 0;
	  pvalue_bonf = 0;
	} else if (pvalue_bonf != pvalue_bonf) {
	  pvalue_bonf = 0;
	}
	    
	search_target_file << scan_num << "\t"; 
	search_target_file << charge << "\t"; 
	search_target_file << precursor_mz << "\t";
	search_target_file << precursor_mass << "\t";
	search_target_file << scores[score_index].second.mass() << "\t";
	search_target_file << scores[score_index].first <<"\t";
	search_target_file << (ntargets+1) << "\t";
	search_target_file << pvalue << "\t";
	search_target_file << target_xpeptides.size() << "\t";
	search_target_file << scores[score_index].second<<endl;
      }
      score_index++;
    }

    delete [] linked_decoy_scores_array;
    //free_spectrum(spectrum);

    cout <<"Done with spectrum"<<endl;
  } // get next spectrum
  search_target_file.close();
  search_decoy_file.close();
  //free_spectrum_collection(spectra);
  //free_spectrum(spectrum);
  return 0;
}

// get all precursor ions within given mass window
void get_ions_from_mz_range(vector<LinkedPeptide>& filtered_ions,
	vector<LinkedPeptide>& all_ions,
	FLOAT_T precursor_mass,
	int charge,
	FLOAT_T mass_window,
	int decoy_iterations) {
  FLOAT_T min_mass = precursor_mass - mass_window;
  FLOAT_T max_mass = precursor_mass + mass_window;
  //cout << "min mass " << min_mass << " max mass " << max_mass << endl;
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


#define BONF_CUTOFF_P 1e-4
#define BONF_CUTOFF_NP 1e-2

double bonf_correct(double nlp_value, int n) {
  if (nlp_value != nlp_value) return 0;
  if (nlp_value == 0) return 0;

  double NL_BONF_CUTOFF_P = (-log(BONF_CUTOFF_P));
  double NL_BONF_CUTOFF_NP= (-log(BONF_CUTOFF_NP));


  double ans = nlp_value - log((double)n);
 
  if ((nlp_value <= NL_BONF_CUTOFF_P) || 
      (ans <= NL_BONF_CUTOFF_NP)) { 
    double p = exp(-nlp_value);
    ans = -log(1-pow((1-p), n));
  }
  return ans;
}
