#include "xhhc_scorer.h"
#include "xhhc_ion_series.h"
//#include "xhhc_search.h"
#include "xlink_compute_qvalues.h"

#include "SearchForXLinks.h"
#include "XLinkMatch.h"
#include "XLinkMatchCollection.h"
#include "XLinkBondMap.h"
#include "XLinkPeptide.h"


//CRUX INCLUDES
#include "objects.h"
#include "FilteredSpectrumChargeIterator.h"
#include "MatchSearch.h"
#include "OutputFiles.h"
#include "SpectrumCollectionFactory.h"


//C++ Includes
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <ctime>



using namespace std;

int SearchForXLinks::xlink_search_main() {

  carp(CARP_INFO, "Beginning crux xlink-search-mods");


  //int num_missed_cleavages = 0;

  /* Get parameters */
  char* ms2_file = get_string_parameter("ms2 file");
  char* input_file = get_string_parameter("protein database");

  /* Prepare input, fasta or index */
  Index* index = NULL;
  Database* database = NULL;
  int num_proteins = prepare_protein_input(input_file, &index, &database);
  carp(CARP_INFO,"Number of proteins:%d",num_proteins);
  free(input_file);

  XLinkBondMap bondmap;
  
  int min_weibull_points = 
    get_int_parameter("min-weibull-points");
  
  int scan_num = 0;

  SpectrumZState zstate;


  //int max_ion_charge = get_max_ion_charge_parameter("max-ion-charge");

  int top_match = get_int_parameter("top-match");
  XLinkPeptide::setLinkerMass(get_double_parameter("link mass"));

  carp(CARP_INFO,"Loading Spectra");
  Spectrum* spectrum = new Spectrum();
  SpectrumCollection* spectra = SpectrumCollectionFactory::create(ms2_file);
  spectra->parse();

  FilteredSpectrumChargeIterator* spectrum_iterator =
    new FilteredSpectrumChargeIterator(spectra);


  /* Prepare output files */
  
  OutputFiles output_files(this);
  output_files.writeHeaders(num_proteins);

  // main loop over spectra in ms2 file
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

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
    XLinkMatchCollection target_candidates(precursor_mz,
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
    if (target_candidates.getMatchTotal() < 1) {
      carp(CARP_INFO, "not enough precursors found in range, skipping scan %d charge %d", scan_num, zstate.getCharge());
      continue;
    }

    //score targets
    carp(CARP_DEBUG, "scoring targets:%d", target_candidates.getMatchTotal());
    target_candidates.scoreSpectrum(spectrum);

    
    carp(CARP_DEBUG,"Getting decoy candidates");

    XLinkMatchCollection decoy_candidates;
    target_candidates.shuffle(decoy_candidates);

    carp(CARP_DEBUG,"scoring decoys");
    decoy_candidates.scoreSpectrum(spectrum);

    if (target_candidates.getMatchTotal() >= min_weibull_points) {
      carp(CARP_INFO,"Fitting weibull to targets");
      target_candidates.fitWeibull();
    //TODO
    //} else if (target_candidates.size() + decoy_candidates.size() >= min_wiebull_points) {
    //  fit_weibull(target_candidates, decoy_candidates, shift, eta, beta, corr);
    } else {
    

      carp(CARP_DEBUG,"Getting weibull training candidates");
      XLinkMatchCollection train_target_candidates(precursor_mz,
                                                   zstate,
                                                   bondmap,
                                                   index,
                                                   database,
                                                   peptide_mods,
                                                   num_peptide_mods,
                                                   TRUE);
   
      XLinkMatchCollection train_candidates(train_target_candidates);
      //get enough weibull training candidates by shuffling.
      carp(CARP_DEBUG,"Shuffling %d:%d", 
        train_candidates.getMatchTotal(), 
        min_weibull_points);
    
      while(train_candidates.getMatchTotal() < min_weibull_points) {
        train_target_candidates.shuffle(train_candidates);
      }
      carp(CARP_DEBUG, "Have %d training candidates", train_candidates.getMatchTotal());
      carp(CARP_DEBUG,"scoring training points");
      train_candidates.scoreSpectrum(spectrum);
      train_candidates.fitWeibull();
    }

    carp(CARP_DEBUG,"setting ranks");
    target_candidates.setRanks();
    decoy_candidates.setRanks();
    
    

    int nprint = min(top_match,target_candidates.getMatchTotal());

    carp(CARP_DEBUG,"Printing %d targets", nprint);
 
    //calculate pvalues.
    for (int idx=0;idx < nprint;idx++) {
      target_candidates.computeWeibullPValue(idx);
    }

    nprint = min(top_match,(int)decoy_candidates.getMatchTotal());

    carp(CARP_DEBUG,"Printing %d decoys", nprint);

    for (int idx=0;idx < nprint;idx++) {
      decoy_candidates.computeWeibullPValue(idx);
    }
    
    vector<MatchCollection*> decoy_vec;
    decoy_vec.push_back(&decoy_candidates);

    //print out
    MatchSearch::printSpectrumMatches(
      output_files,
      (MatchCollection*)&target_candidates, 
      decoy_vec, 
      spectrum,
      false,
      1);
 
    carp(CARP_DEBUG,"Cleanup");
    XLink::deleteAllocatedPeptides();
    
    //carp(CARP_DEBUG, "num targets:%d",target_candidates.size());
    //free_spectrum(spectrum);

    carp(CARP_DEBUG,"Done with spectrum %d", scan_num);
  } // get next spectrum

  output_files.writeFooters();

  // clean up
  delete spectrum_iterator;
  delete spectra;
  for(int mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++){
    free_peptide_mod(peptide_mods[mod_idx]);
  }
  free(peptide_mods);
  Index::free(index);
  Database::freeDatabase(database);

  //Calculate q-values.
  carp(CARP_DEBUG, "Computing Q-Values");
  xlink_compute_qvalues();

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux search-for-xlink-mods.");

  return(0);
}

