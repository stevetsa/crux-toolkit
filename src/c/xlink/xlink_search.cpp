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
#include "FilteredSpectrumChargeIterator.h"


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

  initialize_run(XLINK_SEARCH_MODS_COMMAND, argument_list, num_arguments,
		 option_list, num_options, argc, argv);
  
  carp(CARP_INFO, "Beginning crux xlink-search-mods");


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
  
  unsigned int min_weibull_points = 
    (unsigned int)get_int_parameter("min-weibull-points");
  
  int scan_num = 0;

  SpectrumZState zstate;


  //int max_ion_charge = get_max_ion_charge_parameter("max-ion-charge");

  int top_match = get_int_parameter("top-match");
  XLinkPeptide::setLinkerMass(get_double_parameter("link mass"));

  carp(CARP_INFO,"Loading Spectra");
  Spectrum* spectrum = new Spectrum();
  SpectrumCollection* spectra = new SpectrumCollection(ms2_file);
  spectra->parse();

  FilteredSpectrumChargeIterator* spectrum_iterator =
    new FilteredSpectrumChargeIterator(spectra);

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

  FLOAT_T shift, eta, beta, corr;

  int search_count = 0;
  // for every observed spectrum 
  while (spectrum_iterator->hasNext()) {

    spectrum = spectrum_iterator->next(zstate);
    scan_num = spectrum->getFirstScan();

    carp(CARP_DEBUG,"count %d scan %d charge %d", search_count, scan_num, zstate.getCharge());

    if (search_count % 10 == 0)
      carp(CARP_INFO,"count %d scan %d charge %d", search_count, scan_num, zstate.getCharge());
    search_count++;

    FLOAT_T precursor_mz = spectrum->getPrecursorMz();

    carp(CARP_DEBUG,"Getting targets");  
    MatchCandidateVector target_candidates(precursor_mz,
                                           zstate,
					   bondmap,
					   index,
					   database,
					   peptide_mods,
					   num_peptide_mods,
					   FALSE);

    target_candidates.setScan(spectrum->getFirstScan());
    /*
    if (target_candidates.tooManyCandidates()) {
      carp(CARP_WARNING, "too many potential candidates in range, skipping scan %d charge %d", scan_num, zstate.getCharge());
      continue;
    }
    */
    if (target_candidates.size() < 1) {
      carp(CARP_INFO, "not enough precursors found in range, skipping scan %d charge %d", scan_num, zstate.getCharge());
      continue;
    }

    //score targets
    carp(CARP_DEBUG, "scoring targets");
    target_candidates.scoreSpectrum(spectrum);

    
    carp(CARP_DEBUG,"Getting decoy candidates");

    MatchCandidateVector decoy_candidates;
    target_candidates.shuffle(decoy_candidates);

    carp(CARP_DEBUG,"scoring decoys");
    decoy_candidates.scoreSpectrum(spectrum);

    if (target_candidates.size() >= min_weibull_points) {
      carp(CARP_INFO,"Fitting weibull to targets");
      target_candidates.fitWeibull(shift, eta, beta, corr);
    //TODO
    //} else if (target_candidates.size() + decoy_candidates.size() >= min_wiebull_points) {
    //  fit_weibull(target_candidates, decoy_candidates, shift, eta, beta, corr);
    } else {
    

      carp(CARP_DEBUG,"Getting weibull training candidates");
      MatchCandidateVector train_target_candidates(precursor_mz,
                                                   zstate,
                                                   bondmap,
                                                   index,
                                                   database,
                                                   peptide_mods,
                                                   num_peptide_mods,
                                                   TRUE);
   
      MatchCandidateVector train_candidates(train_target_candidates);
      //get enough weibull training candidates by shuffling.
      carp(CARP_DEBUG,"Shuffling %d:%d", train_candidates.size(), min_weibull_points);
    
      while(train_candidates.size() < min_weibull_points) {
        train_target_candidates.shuffle(train_candidates);
      }
      carp(CARP_DEBUG, "Have %d training candidates", train_candidates.size());
      carp(CARP_DEBUG,"scoring training points");
      train_candidates.scoreSpectrum(spectrum);
      train_candidates.fitWeibull(shift, eta, beta, corr);
    }

    carp(CARP_DEBUG,"setting ranks");
    target_candidates.setRanks();
    decoy_candidates.setRanks();
    
    

    int nprint = min(top_match,(int)target_candidates.size());

    carp(CARP_DEBUG,"Printing %d targets", nprint);
 
    //print out data.
    for (int idx=0;idx < nprint;idx++) {
      target_candidates[idx]->computeWeibullPvalue(shift, eta, beta);
      target_file << target_candidates[idx]->getResultString() << endl;
    }

    nprint = min(top_match,(int)decoy_candidates.size());

    carp(CARP_DEBUG,"Printing %d decoys", nprint);

    for (int idx=0;idx < nprint;idx++) {
      decoy_candidates[idx]->computeWeibullPvalue(shift, eta, beta);
      decoy_file << decoy_candidates[idx]->getResultString() << endl;
    }
    
    carp(CARP_DEBUG,"Cleanup");
    XLink::deleteAllocatedPeptides();
    
    //carp(CARP_DEBUG, "num targets:%d",target_candidates.size());
    //free_spectrum(spectrum);

    carp(CARP_DEBUG,"Done with spectrum %d", scan_num);
  } // get next spectrum
  //carp(CARP_INFO,"Done searching spectra");
  target_file.close();
  decoy_file.close();
  delete spectra;
  //free_spectrum(spectrum);
  //carp(CARP_INFO,"Done freeing spectrum collection");

  //Calculate q-values.
  carp(CARP_DEBUG, "Computing Q-Values");
  xlink_compute_qvalues();

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux search-for-xlink-mods.");

  return(0);
}

