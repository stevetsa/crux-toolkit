#include "xhhc_scorer.h"
#include "xhhc_ion_series.h"
//#include "xhhc_search.h"
#include "xlink_compute_qvalues.h"

#include "MatchCandidate.h"
#include "MatchCandidateVector.h"
#include "XLinkBondMap.h"
#include "XLinkPeptide.h"


//CRUX INCLUDES
#include "objects.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <ctime>



using namespace std;

int xlink_search_main(int argc, char** argv) {

  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);

  /* Define optional command line arguments */
  const char* option_list[] = {
    "verbosity",
    "parameter-file",
    "overwrite",
    "output-dir",
    "precursor-window",
    "precursor-window-type",
    "precursor-window-decoy",
    "precursor-window-type-decoy",
    "max-ion-charge",
    "min-weibull-points",
    "xlink-prevents-cleavage",
    "spectrum-min-mass",
    "spectrum-max-mass",
    "spectrum-charge",
    "top-match",
    "xlink-include-linears",
    "xlink-include-deadends",
    "xlink-include-selfloops",
    "xcorr-use-flanks",
    "use-mgf"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {
    "ms2 file", 
    "protein database", 
    "link sites", 
    "link mass"
  };

  int num_arguments = sizeof(argument_list) / sizeof(char*);

  initialize_run(XLINK_SEARCH_COMMAND, argument_list, num_arguments,
		 option_list, num_options, argc, argv);
  
  carp(CARP_INFO, "Beginning crux xlink-search");


  //int num_missed_cleavages = 0;
  char* ms2_file = get_string_parameter("ms2 file");

  //FLOAT_T precursor_window = get_double_parameter("precursor-window");
  //FLOAT_T precursor_window_decoy = get_double_parameter("precursor-window-decoy");
  //WINDOW_TYPE_T precursor_window_type = 
  //  get_window_type_parameter("precursor-window-type");
  //WINDOW_TYPE_T window_type_decoy = 
  //  get_window_type_parameter("precursor-window-type-decoy");

  char* input_file = get_string_parameter("protein database");
  
  INDEX_T* index = NULL;
  DATABASE_T* database = NULL;
  int num_proteins = prepare_protein_input(input_file, &index, &database);
  carp(CARP_INFO,"Number of proteins:%d",num_proteins);
  free(input_file);

  XLinkBondMap bondmap;
  /*
  unsigned int min_weibull_points = 
    (unsigned int)get_int_parameter("min-weibull-points");
  */
  int scan_num = 0;
  int charge = 1;

  //int max_ion_charge = get_max_ion_charge_parameter("max-ion-charge");

  int top_match = get_int_parameter("top-match");
  XLinkPeptide::setLinkerMass(get_double_parameter("link mass"));

  carp(CARP_INFO,"Loading Spectra");
  SPECTRUM_T* spectrum = allocate_spectrum();
  SPECTRUM_COLLECTION_T* spectra = new_spectrum_collection(ms2_file);
  parse_spectrum_collection(spectra);
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* spectrum_iterator = 
	new_filtered_spectrum_charge_iterator(spectra);

  char* output_dir = get_string_parameter("output-dir");

  string target_filename = string(output_dir)+"/search.target.txt";
  string decoy_filename = string(output_dir)+"/search.decoy.txt";
  ofstream target_file(target_filename.c_str());
  ofstream decoy_file(decoy_filename.c_str());

  target_file << MatchCandidate::getResultHeader()<<endl;
  decoy_file  << MatchCandidate::getResultHeader()<<endl;


  // main loop over spectra in ms2 file
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

  int search_count = 0;
  // for every observed spectrum 
  while (filtered_spectrum_charge_iterator_has_next(spectrum_iterator)) {
    //cerr<<"Getting next spectrum"<<endl;
    spectrum = filtered_spectrum_charge_iterator_next(spectrum_iterator, &charge);
    //SCORER_T* scorer = new_scorer(XCORR);
    scan_num = get_spectrum_first_scan(spectrum);

    if (search_count % 1 == 0)
      carp(CARP_INFO,"count %d scan %d charge %d", search_count, scan_num, charge);
    search_count++;

    FLOAT_T precursor_mz = get_spectrum_precursor_mz(spectrum);
    //FLOAT_T precursor_mass = get_spectrum_neutral_mass(spectrum, charge); 
 
    
    MatchCandidateVector target_candidates(precursor_mz, 
					   charge, 
					   bondmap,
					   index,
					   database,
					   peptide_mods,
					   num_peptide_mods,
					   FALSE);

    target_candidates.setScan(get_spectrum_first_scan(spectrum));


    if (target_candidates.size() < 1) {
      carp(CARP_INFO, "not enough precursors found in range, skipping scan %d charge %d", scan_num, charge);
      continue;
    }

    //cerr <<"Found "<<target_candidates.size()<<" targets"<<endl;
    

    MatchCandidateVector decoy_candidates;
    target_candidates.shuffle(decoy_candidates);

    MatchCandidateVector train_target_candidates(precursor_mz,
					         charge,
						 bondmap,
						 index,
						 database,
						 peptide_mods,
						 num_peptide_mods,
						 TRUE);
   
    //cerr<<"Have "<<train_target_candidates.size()<<" target training candidates"<<endl;
    MatchCandidateVector train_candidates(train_target_candidates);
    //get enough weibull training candidates by shuffling.
    //cerr<<"Shuffling "<<min_weibull_points<<endl;
    /*
    while(train_candidates.size() < min_weibull_points) {
      train_target_candidates.shuffle(train_candidates);
    }
    */
    //cerr<<"Have "<<train_candidates.size()<<" training candidates"<<endl;

    //cerr <<"scoring targets"<<endl;
    target_candidates.scoreSpectrum(spectrum);
    //cerr <<"scoring decoys"<<endl;
    decoy_candidates.scoreSpectrum(spectrum);
    //cerr <<"scoring training candidates"<<endl;
    train_candidates.scoreSpectrum(spectrum);

    //cerr <<"Sorting by XCorr"<<endl;
    target_candidates.sortByXCorr();
    decoy_candidates.sortByXCorr();
    
    FLOAT_T shift, eta, beta, corr;

    //cerr<<"Fitting Weibull"<<endl;
    train_candidates.fitWeibull(shift, eta, beta, corr);
   
    int nprint = min(top_match,(int)target_candidates.size());
 
    //print out data.
    for (int idx=0;idx < nprint;idx++) {
      target_candidates[idx]->computeWeibullPvalue(shift, eta, beta);
      target_file << target_candidates[idx]->getResultString() << endl;
    }

    nprint = min(top_match,(int)decoy_candidates.size());

    for (int idx=0;idx < nprint;idx++) {
      decoy_candidates[idx]->computeWeibullPvalue(shift, eta, beta);
      decoy_file << decoy_candidates[idx]->getResultString() << endl;
    }
    
    XLinkPeptide::cleanUp();
    
    //carp(CARP_DEBUG, "num targets:%d",target_candidates.size());
    //free_spectrum(spectrum);

    //carp(CARP_INFO,"Done with spectrum %d", scan_num);
  } // get next spectrum
  //carp(CARP_INFO,"Done searching spectra");
  target_file.close();
  decoy_file.close();
  free_spectrum_collection(spectra);
  //free_spectrum(spectrum);
  //carp(CARP_INFO,"Done freeing spectrum collection");
  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux search-for-xlinks.");

  return(0);
}

