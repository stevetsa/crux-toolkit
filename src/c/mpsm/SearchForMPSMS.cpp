#include "SearchForMPSMS.h"
#include "MPSM_MatchCollection.h"
#include "MPSM_OutputFiles.h"

#include "FilteredSpectrumChargeIterator.h"

#include "SearchProgress.h"

#include "SpectrumCollection.h"

using namespace std;

SearchForMPSMS::SearchForMPSMS() {

}

SearchForMPSMS::~SearchForMPSMS() {
}


int SearchForMPSMS::main(int argc, char** argv) {


    /* Define optional command line arguments */
  const char* option_list[] = {
    "verbosity",
    "version",
    "parameter-file",
    "overwrite",
    "compute-p-values",
    "spectrum-min-mass",
    "spectrum-max-mass",
    "spectrum-charge",
    "scan-number",
    "output-dir",
    "fileroot",
    "num-decoys-per-target",
    "decoy-location",
    "precursor-window",
    "precursor-window-type",
    "mpsm-max-peptides",
    "rtime-threshold",
    "rtime-all2-threshold",
    "rtime-all3-threshold",
    "rtime-default-threshold",
    "mpsm-top-n",
    "mpsm-do-sort",
    "rtime-predictor",
    "mpsm-top-match"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  

  /* Define required command line arguments */
  const char* argument_list[] = {"ms2 file", "protein database"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);


  initialize_run(SEARCH_MPSMS_COMMAND, argument_list, num_arguments,
                 option_list, num_options, argc, argv);

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  // Get input: ms2 file 
  const char* ms2_file = get_string_parameter_pointer("ms2 file");

  // open ms2 file
  SpectrumCollection* spectra = new SpectrumCollection(ms2_file);

  rtime_predictor_ = RetentionPredictor::createRetentionPredictor();

  rtime_threshold_ = get_boolean_parameter("rtime-threshold");
  rtime_all2_threshold_ = get_double_parameter("rtime-all2-threshold");
  rtime_all3_threshold_ = get_double_parameter("rtime-all3-threshold");
  rtime_default_threshold_ = get_double_parameter("rtime-default-threshold");


  // parse the ms2 file for spectra
  carp(CARP_INFO, "Reading in ms2 file %s", ms2_file);
  if(!spectra->parse()){
    carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2_file);
  }
  
  carp(CARP_DEBUG, "There were %i spectra found in the ms2 file",
       spectra->getNumSpectra());

  /* Get input: protein file */
  char* input_file = get_string_parameter("protein database");

  /* Prepare input, fasta or index */
  INDEX_T* index = NULL;
  DATABASE_T* database = NULL;
  int num_proteins = prepare_protein_input(input_file, &index, &database); 
  free(input_file);
  /*
  char* decoy_input_file = "./shuffled_index";
  INDEX_T* index_shuffled = NULL;
  DATABASE_T* database_shuffled = NULL;
  int num_proteins_shuffled = 
    prepare_protein_input(decoy_input_file, &index_shuffled, &database_shuffled);
  */
  carp(CARP_DEBUG, "Found %i proteins", num_proteins);
  if( num_proteins == 0 ){
    carp(CARP_FATAL, "No proteins were found in the protein source.");
  }
  
  /* Prepare output files */
  MPSM_OutputFiles output_files(SEARCH_COMMAND); 
  output_files.writeHeaders(num_proteins);
  // TODO (BF oct-21-09): consider adding pvalue file to OutputFiles
  FILE* decoy_pvalue_file = NULL;
  if( get_boolean_parameter("decoy-p-values") ){
    carp(CARP_DEBUG, "Opening decoy p-value file.");
    char* decoy_pvalue_filename 
      = get_string_parameter("search-decoy-pvalue-file");
    prefix_fileroot_to_name(&decoy_pvalue_filename);
    char* output_directory = get_string_parameter("output-dir");
    decoy_pvalue_file = create_file_in_path(decoy_pvalue_filename, 
                                            output_directory, 
                                            get_boolean_parameter("overwrite"));
    free(decoy_pvalue_filename);
    free(output_directory);
  }

  /* Perform search: loop over spectra*/
  wall_clock();
  // create spectrum iterator
  FilteredSpectrumChargeIterator* spectrum_iterator = 
    new FilteredSpectrumChargeIterator(spectra);

  // get search parameters for match_collection
  BOOLEAN_T compute_pvalues = get_boolean_parameter("compute-p-values");
  BOOLEAN_T combine_target_decoy = get_boolean_parameter("tdc");
  int num_decoy_files = get_int_parameter("num-decoy-files");

  // For remembering and reporting number of searches
  SearchProgress progress;

  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

  Spectrum* current_spectrum = NULL;

  MPSM_ZStateMap spsm_map;
  MPSM_ZStateMap mpsm_map;

  // for each spectrum
  while(spectrum_iterator->hasNext()){
    SpectrumZState zstate;
    Spectrum* spectrum = 
      spectrum_iterator->next(zstate);
    bool is_decoy = false;

    //we want to collect all charges for a particular spectrum.
    if (current_spectrum == NULL) {
      current_spectrum = spectrum;
    }

    carp(CARP_INFO, 
      "processing spec %d charge:%d mass:%f", 
      spectrum->getFirstScan(), 
      zstate.getCharge(),
      zstate.getNeutralMass());

    if (spectrum != current_spectrum) {
      carp(CARP_INFO, 
        "Processed all charges for spec %d", 
        current_spectrum->getFirstScan());
      carp(CARP_INFO, "Searching for mpsms");
      search(spsm_map, mpsm_map);
      mpsm_map.calcXCorrRanks();
      if (get_boolean_parameter("mpsm-do-sort")) {
        //spsm_map.calcDeltaCN();
        //spsm_map.calcZScores();
        //spsm_map.sortMatches(XCORR);
        //cerr<<"Calculating delta cn"<<endl;
        mpsm_map.calcDeltaCN();
        //cerr<<"Calculating zscores"<<endl;
        mpsm_map.calcZScores();
        //cerr<<"Calculating xcorr ranks"<<endl;
        
      }
      //print out map
      //output the spsms.
      //output_files.writeMatches(spsm_map);
      //cerr<< "writing matches";
      carp(CARP_INFO, "writing matches:%d",mpsm_map.size());
      output_files.writeMatches(mpsm_map);
      //cerr<<"Clear map"<<endl;
      //clear map and clean up match collections.
      mpsm_map.clearMap();

      //cerr<<"Deleting spsm matches"<<endl;
      for (MPSM_ZStateMap::iterator iter = spsm_map.begin();
        iter != spsm_map.end();
        ++iter) {

        vector<MPSM_MatchCollection>& spsm_match_collections = iter -> second;
        for (int idx = 0;idx < spsm_match_collections.size();idx++) {
          spsm_match_collections.at(idx).free();
        }
      }

      spsm_map.clearMap();

      //carp(CARP_FATAL,"Stopping here for now");
      current_spectrum = spectrum;
    }
    

    progress.report(spectrum->getFirstScan(), zstate.getCharge());
 
    // with the target database decide how many peptide mods to use
    MATCH_COLLECTION_T* target_psms = new_empty_match_collection(is_decoy); 
    int max_pep_mods = searchPepMods( target_psms, 
                                        is_decoy,   
                                        index,       
                                        database, 
                                        spectrum, 
                                        zstate,
                                        peptide_mods, 
                                        num_peptide_mods,
                                        compute_pvalues); 
 

    // are there any matches?
    if( get_match_collection_match_total(target_psms) == 0 ){
      // don't print and don't search decoys
      carp(CARP_WARNING, "No matches found for spectrum %i, charge %i",
           spectrum->getFirstScan(), zstate.getCharge());
      free_match_collection(target_psms);
      progress.increment(FALSE);
      continue; // next spectrum
    }

    // now search decoys with the same number of mods
    is_decoy = TRUE;
    // create separate decoy match_collections
    int num_decoy_collections = get_int_parameter("num-decoys-per-target"); 
    MATCH_COLLECTION_T** decoy_collection_list = 
      (MATCH_COLLECTION_T**)mycalloc(sizeof(MATCH_COLLECTION_T*), 
                                     num_decoy_collections);

    int decoy_idx = 0;
    for(decoy_idx = 0; decoy_idx < num_decoy_collections; decoy_idx++){

      MATCH_COLLECTION_T* decoy_psms = new_empty_match_collection(is_decoy);
      decoy_collection_list[decoy_idx] = decoy_psms;

      searchPepMods(decoy_psms, 
                      is_decoy,   
                      index, 
                      database, 
                      spectrum, 
                      zstate, 
                      peptide_mods, 
                      max_pep_mods,
                      compute_pvalues);
    }

   
    //add all matches to an mpsm collection.
    
    //add to map, indexed by charge

    
    vector<MPSM_MatchCollection> mpsm_collections;


    MPSM_MatchCollection spsm_targets(target_psms);
    mpsm_collections.push_back(spsm_targets);
    
    for (decoy_idx = 0; decoy_idx < num_decoy_collections; decoy_idx++) {
      MPSM_MatchCollection spsm_decoy(decoy_collection_list[decoy_idx]);
      mpsm_collections.push_back(spsm_decoy);
    }
    
    spsm_map.insert(mpsm_collections);

   

    progress.increment(TRUE);

    // clean up
    
    //free_match_collection(target_psms);
    for(decoy_idx = 0; decoy_idx < num_decoy_collections; decoy_idx++){
      decoy_collection_list[decoy_idx] = NULL;
    }
    //free the array.
    free(decoy_collection_list);

  }// next spectrum

  //process last spectrum.
  search(spsm_map, mpsm_map);
  mpsm_map.calcXCorrRanks();
  if (get_boolean_parameter("mpsm-do-sort")) { 
    mpsm_map.sortMatches(XCORR);
  //print out map
  //output the spsms.
    mpsm_map.calcDeltaCN();
    mpsm_map.calcZScores();
  }
  //output_files.writeMatches(spsm_map);
  output_files.writeMatches(mpsm_map);

  //clear map and clean up match collections.
  mpsm_map.clearMap();

  for (MPSM_ZStateMap::iterator iter = spsm_map.begin();
    iter != spsm_map.end();
    ++iter) {
    
    vector<MPSM_MatchCollection>& spsm_match_collections = iter -> second;
    for (int idx = 0;idx < spsm_match_collections.size();idx++) {
      spsm_match_collections.at(idx).free();
    }
  }

  spsm_map.clearMap();

  // Finished Searching!

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux-search-for-mpsms");
  free(spectrum_iterator);
  delete spectra;
  free_parameters();

  exit(0);
}// end main

string SearchForMPSMS::getName() {
  return "search-for-mpsms";
}

string SearchForMPSMS::getDescription() {
  return 
    "Search a collection of spectra against a sequence "
    "database returning a collection of matches "
    "corresponding to peptide matches allowing more than one peptide to "
    "be assigned to a spectrum "
    "scored by XCorr.";
}




/* Private functions */

void SearchForMPSMS::search(
  MPSM_ZStateMap& charge_spsm_map, 
  MPSM_ZStateMap& charge_mpsm_map
) {

  charge_mpsm_map.clear();
  
  charge_mpsm_map = charge_spsm_map;

  int max_peptides = min(charge_spsm_map.size(), (size_t)get_int_parameter("mpsm-max-peptides"));

  carp(CARP_INFO,"Max possible peptides/spectrum:%d",max_peptides);

  for (int npeptides = 2; npeptides <= max_peptides;npeptides++) {

    carp(CARP_INFO, "Finding all %d-psm matches", npeptides);
    bool added = extendChargeMap(charge_spsm_map,
      charge_mpsm_map,
      npeptides-1);

    if (!added) {
      break;
    }

  }


}

bool SearchForMPSMS::passRTimeThreshold(
  MPSM_Match& match
) {

  FLOAT_T rtime_max_diff = rtime_predictor_ -> calcMaxDiff(match);
  ZStateIndex& zstate_index = match.getZStateIndex();

  
  
  bool ans = false;

  if (rtime_threshold_) {
    FLOAT_T fdiff = fabs(rtime_max_diff);
    if (zstate_index.numCharge(2) == match.numMatches()) {
      ans = fdiff <= rtime_all2_threshold_;
    } else if (zstate_index.numCharge(3) == match.numMatches()) {
      ans = fdiff <= rtime_all3_threshold_;
    } else {
      ans = fdiff <= rtime_default_threshold_;
    }
  } else {
    ans = true;
  }

  if (ans) {
    match.setRTimeMaxDiff(rtime_max_diff);
  }

  return ans;

}

bool SearchForMPSMS::extendMatch(
  MPSM_Match& orig_mpsm,
  MPSM_MatchCollection& spsm_matches,
  MPSM_ZStateMap& new_mpsm_matches,
  int match_collection_idx) {


  bool match_added = false;
  
  for (int idx=0;idx < spsm_matches.numMatches();idx++) {

    MPSM_Match& match_to_add = spsm_matches[idx];

    MPSM_Match new_match(orig_mpsm);
    bool canadd = new_match.addMatch(match_to_add.getMatch(0));
    if (!canadd) {
      continue;
    }

    bool visited = new_mpsm_matches.visited(new_match, match_collection_idx);
    if (visited) {
      continue;
    }

    bool pass_threshold = passRTimeThreshold(new_match);

    if (pass_threshold) {
      //cerr<<new_match<<":"<<new_match.getRTimeMaxDiff()<<endl;
      new_mpsm_matches.insert(new_match, match_collection_idx);
      match_added = true;
    }
  }

  return match_added;

}
  

bool SearchForMPSMS::extendChargeMap(
  MPSM_ZStateMap& spsm_map,
  MPSM_ZStateMap& current_mpsm_map,
  int mpsm_level) {

  bool success = false;

  MPSM_ZStateMap::iterator map_iter;
  MPSM_ZStateMap::iterator map_iter2;

  set<ZStateIndex> new_zstates;

  for (map_iter = current_mpsm_map.begin();
    map_iter != current_mpsm_map.end();
    ++map_iter) {

    ZStateIndex zstate_index = map_iter->first;
    
    carp(CARP_INFO,"mpsm_level:%d zstate size:%d",mpsm_level, zstate_index.size());
    
    if (zstate_index.size() == mpsm_level) {
    
      //cerr << "Extending zstate:"<<zstate_index<<endl;
        
      vector<MPSM_MatchCollection>& mpsm_match_collections = map_iter -> second;
      MPSM_MatchCollection& mpsm_target_match_collection =  mpsm_match_collections.at(0);

      for (map_iter2 = spsm_map.begin();
        map_iter2 != spsm_map.end();
        ++map_iter2) {
        ZStateIndex zstate_index2 = map_iter2->first;

        ZStateIndex new_zstate_index(zstate_index);

        if ((new_zstate_index.add(zstate_index2.at(0))) &&
            (current_mpsm_map.find(new_zstate_index) == current_mpsm_map.end())) {

          //cerr << "With:"<<zstate_index2<<endl;
          //cerr << "New:"<<new_zstate_index<<endl;

          current_mpsm_map.insert(new_zstate_index);

          vector<MPSM_MatchCollection>& spsm_match_collections = map_iter2 -> second;
          MPSM_MatchCollection& spsm_target_match_collection = spsm_match_collections.at(0);

          for (int mpsm_idx = 0;
            mpsm_idx < mpsm_target_match_collection.numMatches();
            mpsm_idx++) {

            MPSM_Match& mpsm_match = mpsm_target_match_collection[mpsm_idx];
            //cerr <<"Extending "<<mpsm_match.getSequenceString() <<endl;
            success |= extendMatch(
              mpsm_match, 
              spsm_target_match_collection,
              current_mpsm_map,
              0);
          }

          if (success) {

            for (int decoy_idx = 1; decoy_idx < spsm_match_collections.size();decoy_idx++) {
              MPSM_MatchCollection& mpsm_decoy_match_collection = mpsm_match_collections.at(decoy_idx);
              MPSM_MatchCollection& spsm_decoy_match_collection = spsm_match_collections.at(decoy_idx);

              for (int mpsm_idx = 0;
                mpsm_idx < mpsm_decoy_match_collection.numMatches();
                mpsm_idx++) {

                MPSM_Match& mpsm_match = mpsm_decoy_match_collection[mpsm_idx];
                success |= extendMatch(
                  mpsm_match,
                  spsm_decoy_match_collection,
                  current_mpsm_map,
                  decoy_idx);
              }

            }

          }



        } /* if ((new_zstate_index */        
      } /* for (map_iter2 */
    } /* if zstate_index.size() == mpsm_level */
  } /* for (map_iter */

  return success;

}


int SearchForMPSMS::searchPepMods(
  MATCH_COLLECTION_T* match_collection, ///< store PSMs here
  BOOLEAN_T is_decoy,   ///< generate decoy peptides from index/db
  INDEX_T* index,       ///< index to use for generating peptides
  DATABASE_T* database, ///< db to use for generating peptides
  Spectrum* spectrum,         ///< spectrum to search
  SpectrumZState& zstate,       ///< seach spectrum at this charge state
  PEPTIDE_MOD_T** peptide_mods, ///< list of peptide mods to apply
  int num_peptide_mods, ///< how many p_mods to use from the list
  BOOLEAN_T store_scores///< keep all scores for p-value estimation
) {
  
  // set match_collection charge
  set_match_collection_zstate(match_collection, zstate);

  // get spectrum precursor mz
  double mz = spectrum->getPrecursorMz();

  int mod_idx = 0;

  // assess scores after all pmods with x amods have been searched
  int cur_aa_mods = 0;

  // for each peptide mod
  for(mod_idx=0; mod_idx<num_peptide_mods; mod_idx++){
    // get peptide mod
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];

    // is it time to assess matches?
    int this_aa_mods = peptide_mod_get_num_aa_mods(peptide_mod);
    
    if( this_aa_mods > cur_aa_mods ){
      carp(CARP_DEBUG, "Finished searching %i mods", cur_aa_mods);
      BOOLEAN_T passes = isSearchComplete(match_collection, cur_aa_mods);
      if( passes ){
        carp(CARP_DETAILED_DEBUG, 
             "Ending search with %i modifications per peptide", cur_aa_mods);
        break;
      }// else, search with more mods
      cur_aa_mods = this_aa_mods;
    }
    
    // get peptide iterator
    MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator =
      new_modified_peptides_iterator_from_zstate(mz,
                                             zstate,
                                             peptide_mod, 
                                             is_decoy,
                                             index,
                                             database);
    
    
    // score peptides
    int added = add_matches(match_collection, 
                            spectrum, 
                            zstate, 
                            peptide_iterator,
                            is_decoy,
                            store_scores,
                            get_boolean_parameter("compute-sp"),
                            FALSE // don't filtery by Sp
                            );
    
    carp(CARP_DEBUG, "Added %i matches", added);

    free_modified_peptides_iterator(peptide_iterator);
    
  }//next peptide mod

  return mod_idx;

}

bool SearchForMPSMS::isSearchComplete(
  MATCH_COLLECTION_T* matches, 
  int mods_per_peptide) {

  if( matches == NULL ){
    return false;
  }

  // keep searching if no limits on how many mods per peptide
  if( get_int_parameter("max-mods") == MAX_PEPTIDE_LENGTH ){
    return false;
  }
  // stop searching if at max mods per peptide
  if( mods_per_peptide == get_int_parameter("max-mods") ){ 
    return true;
  }

  // test for minimun score found

  return false;

}

