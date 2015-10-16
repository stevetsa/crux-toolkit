/**
 * \file xlink_search.cpp
 * \brief Object for running search-for-xlinks (new code)
 *****************************************************************************/
#include "SearchForXLinks.h"
#include "XLinkMatch.h"
#include "XLinkMatchCollection.h"
#include "XLinkBondMap.h"
#include "XLinkPeptide.h"
#include "XLinkIonSeriesCache.h"
#include "xlink_compute_qvalues.h"

//CRUX INCLUDES
#include "objects.h"
#include "model/FilteredSpectrumChargeIterator.h"
#include "io/OutputFiles.h"
#include "io/SpectrumCollectionFactory.h"
#include "XLinkDatabase.h"


//C++ Includes
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <ctime>



using namespace std;

void buildArguments(
  vector<string>& args_vec, 
  int &argc, 
  char** &argv) {

  argc = args_vec.size();
  argv = new char*[argc];

  argv[0] = (char*)args_vec[0].c_str();
  carp(CARP_INFO, "argv[0]=%s", argv[0]);
  for (int idx = 1; idx < argc; idx++) {
    argv[idx] = (char*)args_vec[idx].c_str();
    carp(CARP_INFO, "argv[%i]=%s", idx, argv[idx]);
  }

}

void writeTrainingCandidates(XLinkMatchCollection* training_candidates, int scan_num) {

  training_candidates->sort(XCORR);

  string output_dir = get_string_parameter("output-dir");

  string output_file = output_dir + "/" +
    "scan." + DelimitedFileWriter::to_string(scan_num) + 
    ".charge." + DelimitedFileWriter::to_string(training_candidates->getCharge()) +
    ".training.candidates.txt";
  cerr<<"writing "<<output_file<<endl;
  DelimitedFileWriter writer(output_file.c_str());

  writer.setColumnName("scan", 0);
  writer.setColumnName("charge", 1);
  writer.setColumnName("decoy", 2);
  writer.setColumnName("sequence", 3);
  writer.setColumnName("eta", 4);
  writer.setColumnName("beta", 5);
  writer.setColumnName("shift", 6);
  writer.setColumnName("xcorr score", 7);
  writer.setColumnName("p-value", 8);
  
  writer.writeHeader();

  for (size_t idx = 0 ;idx < training_candidates->getMatchTotal();idx++) {
    writer.setColumnCurrentRow(0, scan_num);
    writer.setColumnCurrentRow(2, (*training_candidates)[idx]->isDecoy());
    writer.setColumnCurrentRow(1, (*training_candidates)[idx]->getCharge());
    writer.setColumnCurrentRow(3, (*training_candidates)[idx]->getSequenceString());
    writer.setColumnCurrentRow(4, training_candidates->getEta());
    writer.setColumnCurrentRow(5, training_candidates->getBeta());
    writer.setColumnCurrentRow(6, training_candidates->getShift());
    writer.setColumnCurrentRow(7, (*training_candidates)[idx]->getScore(XCORR));
    training_candidates->computeWeibullPValue(idx);
    writer.setColumnCurrentRow(8, (*training_candidates)[idx]->getPValue());
    

    writer.writeRow();
  }

}


/**
 * main method for SearchForXLinks that implements to refactored code
 */
int SearchForXLinks::xlinkSearchMain() {
  
  carp(CARP_INFO, "Beginning crux search-for-xlinks (new code)");

  /* Get parameters */
  carp(CARP_INFO, "Getting parameters");
  string ms2_file = get_string_parameter("ms2 file");
  string input_file = get_string_parameter("protein fasta file");
  string output_directory = get_string_parameter("output-dir");
  int top_match = get_int_parameter("top-match");
  XLinkPeptide::setLinkerMass(get_double_parameter("link mass"));
  int min_weibull_points = get_int_parameter("min-weibull-points");
  bool compute_pvalues = get_boolean_parameter("compute-p-values");

  XLinkBondMap bondmap;

  /* Prepare input fasta  */


  carp(CARP_DEBUG, "Preparing database");
  XLinkDatabase::initialize();
  Database* database = NULL;
  int num_proteins = 0;//prepare_protein_input(input_file, &database);
  //carp(CARP_DEBUG, "Number of proteins:%d",num_proteins);
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = 0;

  /* Usually for debugging purposes, print out the database of canddiates */
  if (get_boolean_parameter("xlink-print-db"))
  {

    carp(CARP_INFO, "generating and printing xlink database");
    XLinkDatabase::print();
    return 0;
  }

  int scan_num = 0;
  int skipped_no_candidates = 0;

  SpectrumZState zstate;

  carp(CARP_DEBUG, "Loading Spectra");
  Crux::Spectrum* spectrum = NULL;
  Crux::SpectrumCollection* spectra = SpectrumCollectionFactory::create(ms2_file);
  spectra->parse();

  FilteredSpectrumChargeIterator* spectrum_iterator =
    new FilteredSpectrumChargeIterator(spectra);

  /* Prepare output files */
  carp(CARP_DEBUG, "Preparing output files");
  OutputFiles output_files(this);
  output_files.writeHeaders(num_proteins);

  // main loop over spectra in ms2 file

  int search_count = 0;
  FLOAT_T num_spectra = (FLOAT_T)spectra->getNumSpectra();

  // for every observed spectrum 
  carp(CARP_DEBUG, "Searching Spectra");


  while (spectrum_iterator->hasNext()) {

    spectrum = spectrum_iterator->next(zstate);
    scan_num = spectrum->getFirstScan();

    carp(CARP_DEBUG,"count %d scan %d charge %d", search_count, scan_num, zstate.getCharge());
    if (search_count > 0 && search_count % 1000 == 0) {
      carp(CARP_INFO, 
	   "%d spectrum-charge combinations searched, %.0f%% complete",
	   search_count, search_count / num_spectra * 100);
    }
    search_count++;

    FLOAT_T precursor_mz = spectrum->getPrecursorMz();

    XLinkMatchCollection* target_candidates = new XLinkMatchCollection(
		       spectrum,
                       zstate,
                       false,
                       false
);

    carp(CARP_DEBUG, "Scan=%d charge=%d mass=%lg candidates=%d", 
      scan_num, 
      zstate.getCharge(), 
      zstate.getNeutralMass(), 
      target_candidates->getMatchTotal());   

    if (target_candidates->getMatchTotal() < 0) {
      carp(CARP_ERROR, "Scan %d has %d candidates.", scan_num, 
	   target_candidates->getMatchTotal());
    } else if (target_candidates->getMatchTotal() == 0) {
      skipped_no_candidates++;
      carp(CARP_DEBUG, "Skipping scan %d charge %d mass %lg", 
        scan_num, 
	zstate.getCharge(),
	zstate.getNeutralMass()
      );
      delete target_candidates;
      //XLink::deleteAllocatedPeptides();
      continue;
    }

    //score targets
    target_candidates->scoreSpectrum(spectrum);
    
    carp(CARP_DEBUG, "Getting decoy candidates");

    //XLinkMatchCollection* decoy_candidates = new XLinkMatchCollection();
    //target_candidates->shuffle(*decoy_candidates);
    
    XLinkMatchCollection* decoy_candidates = new XLinkMatchCollection(
								      spectrum,
								      zstate,
								      true,
								      false);

    
    carp(CARP_DEBUG, "scoring decoys");
    decoy_candidates->scoreSpectrum(spectrum);

    if (compute_pvalues) {

      XLinkMatchCollection* weibull_set = NULL;

      if (target_candidates->getMatchTotal() >= min_weibull_points) {
        carp(CARP_DEBUG, "Fitting weibull to targets");
        target_candidates->fitWeibull();
        if (get_boolean_parameter("write-weibull-points")) {
	  writeTrainingCandidates(target_candidates, scan_num);
	}
        MatchCollection::transferWeibull(target_candidates, decoy_candidates);
      //TODO
      //} else if (target_candidates.size() + decoy_candidates.size() >= min_wiebull_points) {
      //  fit_weibull(target_candidates, decoy_candidates, shift, eta, beta, corr);
      } else {
    
        carp(CARP_DEBUG, "Getting weibull training candidates");
        XLinkMatchCollection* train_target_candidates =
          new XLinkMatchCollection(
	    spectrum,
            zstate,
            false,
            true);
        
	if (train_target_candidates->getMatchTotal() == 0) {
	  carp(CARP_WARNING, "No candidates found in decoy window?  Getting target");
	  delete train_target_candidates;
	  train_target_candidates = new XLinkMatchCollection(*target_candidates);
	}
        XLinkMatchCollection* train_candidates = new XLinkMatchCollection(*train_target_candidates);
        //get enough weibull training candidates by shuffling.
        carp(CARP_DEBUG, "Shuffling %d:%d", 
          train_candidates->getMatchTotal(), 
          min_weibull_points);
    
        while(train_candidates->getMatchTotal() < min_weibull_points) {
          train_target_candidates->shuffle(*train_candidates);
        }
        carp(CARP_DEBUG, "Have %d training candidates", train_candidates->getMatchTotal());
        carp(CARP_DEBUG, "scoring training points");
        train_candidates->scoreSpectrum(spectrum);
        carp(CARP_DEBUG, "fitting weibul to training points");
        train_candidates->fitWeibull();
        carp(CARP_DEBUG, "transferring weibull parameters");
        MatchCollection::transferWeibull(train_candidates, target_candidates);
        MatchCollection::transferWeibull(train_candidates, decoy_candidates);
        carp(CARP_DEBUG, "cleanup train");
	if (get_boolean_parameter("write-weibull-points")) {
	  writeTrainingCandidates(train_candidates, scan_num);
	}

        delete train_candidates;
        carp(CARP_DEBUG, "cleanup train_target");
        delete train_target_candidates;
      }

      target_candidates->sort(XCORR);


      int nprint = min(top_match,target_candidates->getMatchTotal());

      carp(CARP_DEBUG, "Calculating %d target p-values", nprint);
 
      //calculate pvalues.
      for (int idx=0;idx < nprint;idx++) {
        target_candidates->computeWeibullPValue(idx);
      }

      nprint = min(top_match,(int)decoy_candidates->getMatchTotal());

      carp(CARP_DEBUG,"Calculating %d decoy p-values", nprint);

      decoy_candidates->sort(XCORR);

      for (int idx=0;idx < nprint;idx++) {
        decoy_candidates->computeWeibullPValue(idx);
      }

    } // if (compute_p_values)

    //print out

    vector<MatchCollection*> decoy_vec;
    decoy_vec.push_back(decoy_candidates);

    if (decoy_candidates->getScoredType(SP) == true) {
      decoy_candidates->populateMatchRank(SP);
    }
    decoy_candidates->populateMatchRank(XCORR);
    decoy_candidates->sort(XCORR);


    carp(CARP_DEBUG, "Ranking");

    if (target_candidates->getScoredType(SP) == true) {
      target_candidates->populateMatchRank(SP);
    }
    target_candidates->populateMatchRank(XCORR);
    target_candidates->sort(XCORR);

    
    carp(CARP_DEBUG, "Writing results");
    output_files.writeMatches(
      (MatchCollection*)target_candidates, 
      decoy_vec,
      XCORR,
      spectrum);

    /* Clean up */
    delete decoy_candidates;
    delete target_candidates;
    //XLink::deleteAllocatedPeptides();
    
    //free_spectrum(spectrum);

    carp(CARP_DEBUG, "Done with spectrum %d", scan_num);
    carp(CARP_DEBUG, "=====================================");
  } // get next spectrum

  carp(CARP_INFO, "Skipped %d (%g%%) spectra with 0 candidates.", 
       skipped_no_candidates, skipped_no_candidates / num_spectra * 100);

  carp(CARP_INFO, "Skipped %d (%g%%) spectra with too few peaks.", 
       spectrum_iterator->numSkipped(), 
       spectrum_iterator->numSkipped() / num_spectra * 100);

  output_files.writeFooters();

  // clean up


  delete spectrum_iterator;
  delete spectra;
  XLink::deleteAllocatedPeptides();
  for(int mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++){
    free_peptide_mod(peptide_mods[mod_idx]);
  }
  free(peptide_mods);

  Scorer::finalize();
  XLinkIonSeriesCache::finalize();
  IonSeries::finalize();
  XLinkDatabase::finalize();
  //Calculate q-values via p-values from weibull fit.
  if (compute_pvalues) {
    carp(CARP_DEBUG, "Computing Q-Values using P-values");
    xlink_compute_qvalues();
  }

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux search-for-xlinks (new).");

  return(0);
}

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
