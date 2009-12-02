#include "xhhc_scorer.h"
#include "xhhc_ion_series.h"
//#include "xhhc_search.h"

//CRUX INCLUDES
extern "C" {
#include "objects.h"
}


#include <fstream>
#include <math.h>

#include <ctime>

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

#define NUM_XLINK_SEARCH_OPTIONS 14
#define NUM_XLINK_SEARCH_ARGS 4

int xlink_search_main(int argc, char** argv) {

  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);

  /* Define optional command line arguments */
  int num_options = NUM_XLINK_SEARCH_OPTIONS;
  const char* option_list[NUM_XLINK_SEARCH_OPTIONS] = {
    "verbosity",
    "version",
    "parameter-file",
    "overwrite",
    "output-dir",
    "mass-window",
    "mass-window-decoy",
    "min-weibull-points",
    "missed-link-cleavage",
    "top-match",
    "xlink-include-linears",
    "xlink-include-deadends",
    "xlink-include-selfloops",
    "xcorr-use-flanks"
  };

  /* Define required command line arguments */
  int num_arguments = NUM_XLINK_SEARCH_ARGS;
  const char* argument_list[NUM_XLINK_SEARCH_ARGS] = {"ms2 file", 
						      "protein input", 
						      "link sites", 
						      "link mass"};

  /* Initialize parameter.c and set default values*/
  initialize_parameters();

  /* Define optional and required arguments */
  select_cmd_line_options(option_list, num_options);
  select_cmd_line_arguments(argument_list, num_arguments);

  /* Parse the command line, including optional params file
     Includes syntax, type, and bounds checking, dies on error */
  parse_cmd_line_into_params_hash(argc, argv, "crux xlink-search");

  /* Set seed for random number generation */
  if(strcmp(get_string_parameter_pointer("seed"), "time")== 0){
    time_t seconds; // use current time to seed
    time(&seconds); // Get value from sys clock and set seconds variable.
    srand((unsigned int) seconds); // Convert seconds to a unsigned int
  }
  else{
    srand((unsigned int)atoi(get_string_parameter_pointer("seed")));
  }
  
  /* Create output directory */ 
  char* output_directory = get_string_parameter("output-dir");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  int result = create_output_directory(
    output_directory, 
    overwrite
  );
  if( result == -1 ){
    carp(CARP_FATAL, "Unable to create output directory %s.", output_directory);
  }

  /* Open the log file to record carp messages */
  char* log_file_name = get_string_parameter("search-log-file");
  open_log_file(&log_file_name);
  free(log_file_name);
  log_command_line(argc, argv);

  // Write the parameter file
  char* param_file_name = get_string_parameter("search-param-file");
  print_parameter_file(&param_file_name);
  free(param_file_name);
  
  carp(CARP_INFO, "Beginning crux xlink-search");

  char* missed_link_cleavage = get_string_parameter("missed-link-cleavage");

  int num_missed_cleavages = 0;
  char* ms2_file = get_string_parameter("ms2 file");

  FLOAT_T mass_window = get_double_parameter("mass-window");
  FLOAT_T mass_window_decoy = get_double_parameter("mass-window-decoy");

  char* database = get_string_parameter("protein input");
  char* links = get_string_parameter("link sites");

  int min_weibull_points = get_int_parameter("min-weibull-points");

  int scan_num = 0;
  int charge = 1;

  int top_match = get_int_parameter("top-match");

  FLOAT_T linker_mass = get_double_parameter("link mass");

  LinkedPeptide::linker_mass = linker_mass;
  vector<LinkedPeptide> all_ions;
  carp(CARP_INFO,"Calling find all precursor ions");
  find_all_precursor_ions(all_ions, links, missed_link_cleavage, database,1);
  carp(CARP_INFO,"Sort");
  // sort filtered ions and decoy ions by mass
  //sort(all_ions.begin(), all_ions.end());

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
  search_target_file << "peptide mass mono\t";
  search_target_file << "peptide mass average\t";
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
  search_decoy_file << "peptide mass mono\t";
  search_decoy_file << "peptide mass average\t";
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

    carp(CARP_INFO,"scan %d", scan_num);
    
    //vector<pair<FLOAT_T, LinkedPeptide> > linked_scores;
    //vector<pair<FLOAT_T, LinkedPeptide> > single_scores;
    vector<pair<FLOAT_T, LinkedPeptide> > scores;

    vector<LinkedPeptide> target_xpeptides;
    vector<LinkedPeptide> target_decoy_xpeptides;
    vector<LinkedPeptide> decoy_train_xpeptides;
    vector<LinkedPeptide> decoy_xpeptides;

    FLOAT_T precursor_mz = get_spectrum_precursor_mz(spectrum);
    FLOAT_T precursor_mass = get_spectrum_neutral_mass(spectrum, charge); 



    carp(CARP_INFO, "finding target xpeptides in mass window...%g", mass_window);
    get_ions_from_mz_range(
	target_xpeptides, // stored in this vector
	all_ions,
	precursor_mass,
	charge,
	mass_window,
	0);

    if (target_xpeptides.size() < 1) {
      carp(CARP_INFO, "not enough precursors found in range, skipping %d", scan_num);
      continue;
    }
    

    carp(CARP_INFO, "finding training xpeptides in decoy mass window..%g", mass_window_decoy);
    get_ions_from_mz_range(
	target_decoy_xpeptides,
	all_ions,
	precursor_mass,
	charge,
	mass_window_decoy,
	0);
    
    carp(CARP_DETAILED_DEBUG, "Creating decoys for target window");
    //create the decoys from the target found in the target_mass_window.
    for (vector<LinkedPeptide>::iterator ion = target_xpeptides.begin();
	 ion != target_xpeptides.end(); ++ion) {
        add_decoys(decoy_xpeptides, *ion);
    }
    
    
    carp(CARP_DETAILED_DEBUG, "Creating decoys for decoy mass window");
    //create the decoys from the target found in the decoy_mass_window.
    while (decoy_train_xpeptides.size() < min_weibull_points) {
      for (vector<LinkedPeptide>::iterator ion = target_decoy_xpeptides.begin();
	   ion != target_decoy_xpeptides.end(); ++ion) {
	add_decoys(decoy_train_xpeptides, *ion);
      }
    }    

    carp(CARP_INFO, "num targets:%d",target_xpeptides.size());
    carp(CARP_INFO, "num decoys:%d", decoy_xpeptides.size());
    carp(CARP_INFO, "num training decoys:%d", decoy_train_xpeptides.size());

    clock_t start_clock = clock();


    // for every ion in the mass window
    carp(CARP_INFO, "Scoring targets");
    for (int i=0;i<target_xpeptides.size();i++) {
      LinkedIonSeries ion_series = LinkedIonSeries(links, charge);
      ion_series.add_linked_ions(target_xpeptides[i]);
      score = hhc_scorer.score_spectrum_vs_series(spectrum, ion_series);
      scores.push_back(make_pair(score, target_xpeptides[i]));
    }

    clock_t target_clock = clock();

    carp(CARP_INFO, "Scoring decoys.");
    for (int i=0;i<decoy_xpeptides.size();i++) {
      LinkedIonSeries ion_series = LinkedIonSeries(links, charge);
      ion_series.add_linked_ions(decoy_xpeptides[i]);
      score = hhc_scorer.score_spectrum_vs_series(spectrum, ion_series);
      scores.push_back(make_pair(score, decoy_xpeptides[i]));
    }


    //use the decoy scores to build the estimator.
    // create arrays to pass to crux's weibull methods
    FLOAT_T* linked_decoy_scores_array = new FLOAT_T[decoy_train_xpeptides.size()+target_xpeptides.size()];

    clock_t decoy_clock = clock();
    carp(CARP_INFO, "scoring training decoys...");
    // score all training decoys
    for (int i=0;i<decoy_train_xpeptides.size();i++) {
      LinkedIonSeries ion_series = LinkedIonSeries(links, charge);
      //ion_series.clear();
      ion_series.add_linked_ions(decoy_train_xpeptides[i]);
      score = hhc_scorer.score_spectrum_vs_series(spectrum, ion_series);
      linked_decoy_scores_array[i] = score;
    }
  
    



    clock_t train_decoy_clock = clock();




    
    for (int i=0;i<scores.size();i++) {
      if (!scores[i].second.is_decoy())
	linked_decoy_scores_array[i+decoy_train_xpeptides.size()] = scores[i].first;
    }
    
    // sort scores
    sort(scores.begin(), scores.end(), greater<pair<FLOAT_T, LinkedPeptide> >());

    clock_t create_array_clock = clock();


   // weibull parameters for candidates
    FLOAT_T eta_linked = 0.0;
    FLOAT_T beta_linked  = 0.0;
    FLOAT_T shift_linked  = 0.0;
    FLOAT_T correlation_linked  = 0.0;

    // fit weibull to decoys

    hhc_estimate_weibull_parameters_from_xcorrs(linked_decoy_scores_array, 
						decoy_train_xpeptides.size(), 
						&eta_linked, &beta_linked, 
						&shift_linked, &correlation_linked, 
						spectrum, charge);
    
    clock_t weibull_clock = clock();
    
    double target_time = (double(target_clock) - double(start_clock)) / CLOCKS_PER_SEC;
    double decoy_time = (double(decoy_clock) - double(target_clock)) / CLOCKS_PER_SEC;
    double train_decoy_time = (double(train_decoy_clock) - double(decoy_clock)) / CLOCKS_PER_SEC;
    double create_array_time =(double(create_array_clock) - double(train_decoy_clock)) / CLOCKS_PER_SEC;
    double weibull_time = (double(weibull_clock) - double(create_array_clock)) / CLOCKS_PER_SEC;


    carp(CARP_INFO,"target:%g",target_time);
    carp(CARP_INFO,"decoy:%g",decoy_time);
    carp(CARP_INFO,"train decoy:%g",train_decoy_time);
    carp(CARP_INFO,"create array:%g",create_array_time);
    carp(CARP_INFO,"weibull:%g",weibull_time);
    


    int ndecoys = 0;
    int ntargets = 0;
    int score_index = 0;

    while (score_index < scores.size() && (ndecoys < top_match || ntargets < top_match)) {
      if (scores[score_index].second.is_decoy() && ndecoys < top_match) {
	double pvalue = compute_weibull_pvalue(scores[score_index].first, eta_linked, beta_linked, shift_linked);
	double pvalue_bonf = pvalue;//bonf_correct(pvalue, decoy_xpeptides.size());
	
	if (pvalue != pvalue) {
	  pvalue = 1;
	  pvalue_bonf = 1;
	} else if (pvalue_bonf != pvalue_bonf) {
	  pvalue_bonf = 1;
	}
	
	search_decoy_file << scan_num << "\t"; 
	search_decoy_file << charge << "\t"; 
	search_decoy_file << precursor_mz << "\t";
	search_decoy_file << precursor_mass << "\t";
	search_decoy_file << scores[score_index].second.mass(MONO) << "\t";
	search_decoy_file << scores[score_index].second.mass(AVERAGE) << "\t";
	search_decoy_file << scores[score_index].first <<"\t";
	search_decoy_file << ndecoys << "\t";
	search_decoy_file << pvalue << "\t";
	search_decoy_file << decoy_xpeptides.size() << "\t";
	search_decoy_file << scores[score_index].second<<endl;

	ndecoys++;
      } else if (ntargets < top_match) {
	ntargets++;

	double pvalue = compute_weibull_pvalue(scores[score_index].first, eta_linked, beta_linked, shift_linked);
	double pvalue_bonf = pvalue;//bonf_correct(pvalue, target_xpeptides.size());
	
	if (pvalue != pvalue) {
	  pvalue = 1;
	  pvalue_bonf = 1;
	} else if (pvalue_bonf != pvalue_bonf) {
	  pvalue_bonf = 1;
	}
	    
	search_target_file << scan_num << "\t"; 
	search_target_file << charge << "\t"; 
	search_target_file << precursor_mz << "\t";
	search_target_file << precursor_mass << "\t";
	search_target_file << scores[score_index].second.mass(MONO) << "\t";
	search_target_file << scores[score_index].second.mass(AVERAGE) << "\t";
	search_target_file << scores[score_index].first <<"\t";
	search_target_file << ntargets << "\t";
	search_target_file << pvalue << "\t";
	search_target_file << target_xpeptides.size() << "\t";
	search_target_file << scores[score_index].second<<endl;
      }
      score_index++;
    }

    delete [] linked_decoy_scores_array;
    //free_spectrum(spectrum);

    carp(CARP_DETAILED_DEBUG,"Done with spectrum %d", scan_num);
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
  carp(CARP_DETAILED_DEBUG,"get_ions_from_mz_range()");
  carp(CARP_DETAILED_DEBUG,"min_mass %g max_mass %g", min_mass, max_mass);

  FLOAT_T ion_mass;
  for (vector<LinkedPeptide>::iterator ion = all_ions.begin();
	ion != all_ions.end(); ++ion) {
    ion->set_charge(charge);
    ion->calculate_mass(get_mass_type_parameter("isotopic-mass"));
    ion_mass = ion->mass(get_mass_type_parameter("isotopic-mass"));
    if (ion_mass >= min_mass && ion_mass <= max_mass) {
      filtered_ions.push_back(*ion);
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
