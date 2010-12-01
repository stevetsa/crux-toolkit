/**
 * \file mpsm_search.cpp
 * BASED ON: original_match_search.c
 * DATE: Aug 19, 2008
 * AUTHOR: Sean McIlwain
 * DESCRIPTION: Main file for crux-search-for-matches.  Given an ms2
 * file and a fasta file or index, compare all spectra to peptides in
 * the fasta file/index and return high scoring matches.  Peptides are
 * determined by parameters for length, mass, mass tolerance, cleavages,
 * modifications. Score first by a preliminary method, keep only the
 * top ranking matches, score those with a second method and re-rank
 * by second score.  Output in binary csm file format or text sqt file
 * format. 
 */
/*
 * Here is the outline for how the new search should work

   for each spectrum
     for each charge state
      for each peptide modification
        create a peptide iterator
        for each peptide
         score peptide/spectra
      if passes criteria, print results and move on
      else next peptide modification  
 */
#include "match_collection.h"
#include "carp.h"
#include "crux-utils.h"
#include "parameter.h"
#include "spectrum_collection.h"
#include <errno.h>
#include "output-files.h"
#include "search-progress.h"

#include "MPSM_Match.h"
#include "MPSM_MatchCollection.h"
#include "MPSM_OutputFiles.h"
#include "MPSM_ChargeMap.h"

#include "RetentionPredictor.h"

/* Private functions */

void search_for_mpsms(
  MPSM_ChargeMap& charge_spsm_map, 
  MPSM_ChargeMap& charge_mpsm_map
);


int mpsm_search_pep_mods(
  MATCH_COLLECTION_T* match_collection, ///< store PSMs here
  BOOLEAN_T is_decoy,   ///< generate decoy peptides from index/db
  INDEX_T* index,       ///< index to use for generating peptides
  DATABASE_T* database, ///< db to use for generating peptides
  SPECTRUM_T* spectrum, ///< spectrum to search
  int charge,           ///< seach spectrum at this charge state
  PEPTIDE_MOD_T** pep_mod_list, ///< list of peptide mods to apply
  int num_peptide_mods, ///< how many p_mods to use from the list
  BOOLEAN_T store_scores///< keep all scores for p-value estimation
);
/*
void add_decoy_scores(
  MATCH_COLLECTION_T* target_psms, ///< add scores to these matches
  SPECTRUM_T* spectrum, ///<
  int charge, ///< 
  INDEX_T* index, ///< search this index if not null
  DATABASE_T* database, ///< search this database if not null
  PEPTIDE_MOD_T** peptitde_mods, ///< list of peptide mods to search
  int num_peptide_mods ///< number of mods in the above array
);
*/
BOOLEAN_T mpsm_is_search_complete(MATCH_COLLECTION_T* matches, 
                             int mods_per_peptide);
void print_spectrum_matches(
  OutputFiles& output_files,       
  MATCH_COLLECTION_T* target_psms, 
  MATCH_COLLECTION_T** decoy_psms,
  int num_decoy_collections,
  SPECTRUM_T* spectrum,             
  BOOLEAN_T combine_target_decoy,
  int num_decoy_files
                   );

#ifdef SEARCH_ENABLED // Discard this code in open source release

BOOLEAN_T rtime_threshold = FALSE;
double rtime_all2_threshold = 8.4;
double rtime_all3_threshold = 14.1;
double rtime_default_threshold = 14.7;

int mpsm_search_main(int argc, char** argv){

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
    "top-match"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  

  /* Define required command line arguments */
  const char* argument_list[] = {"ms2 file", "protein database"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  initialize_run(SEARCH_MPSMS_COMMAND, argument_list, num_arguments,
                 option_list, num_options, argc, argv);

  // Get input: ms2 file 
  const char* ms2_file = get_string_parameter_pointer("ms2 file");

  // open ms2 file
  SPECTRUM_COLLECTION_T* spectra = new_spectrum_collection(ms2_file);

  rtime_threshold = get_boolean_parameter("rtime-threshold");
  rtime_all2_threshold = get_double_parameter("rtime-all2-threshold");
  rtime_all3_threshold = get_double_parameter("rtime-all3-threshold");
  rtime_default_threshold = get_double_parameter("rtime-default-threshold");


  // parse the ms2 file for spectra
  carp(CARP_INFO, "Reading in ms2 file %s", ms2_file);
  if(!parse_spectrum_collection(spectra)){
    carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2_file);
  }
  
  carp(CARP_DEBUG, "There were %i spectra found in the ms2 file",
       get_spectrum_collection_num_spectra(spectra));

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
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* spectrum_iterator = 
    new_filtered_spectrum_charge_iterator(spectra);

  // get search parameters for match_collection
  BOOLEAN_T compute_pvalues = get_boolean_parameter("compute-p-values");
  BOOLEAN_T combine_target_decoy = get_boolean_parameter("tdc");
  int num_decoy_files = get_int_parameter("num-decoy-files");

  // For remembering and reporting number of searches
  SearchProgress progress;

  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );

  SPECTRUM_T* current_spectrum = NULL;

  MPSM_ChargeMap spsm_map;
  MPSM_ChargeMap mpsm_map;

  // for each spectrum
  while(filtered_spectrum_charge_iterator_has_next(spectrum_iterator)){
    int charge = 0;
    SPECTRUM_T* spectrum = 
      filtered_spectrum_charge_iterator_next(spectrum_iterator, &charge);
    BOOLEAN_T is_decoy = FALSE;

    //we want to collect all charges for a particular spectrum.
    if (current_spectrum == NULL) {
      current_spectrum = spectrum;
    }

    carp(CARP_DEBUG, "processing spec %d charge:%d", get_spectrum_first_scan(spectrum), charge);
    if (spectrum != current_spectrum) {
      carp(CARP_DEBUG, "Processed all charges for spec %d", get_spectrum_first_scan(current_spectrum));
      carp(CARP_DEBUG, "Searching for mpsms");
      search_for_mpsms(spsm_map, mpsm_map);
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
      output_files.writeMatches(mpsm_map);
      //cerr<<"Clear map"<<endl;
      //clear map and clean up match collections.
      mpsm_map.clearMap();

      //cerr<<"Deleting spsm matches"<<endl;
      for (MPSM_ChargeMap::iterator iter = spsm_map.begin();
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
    

    progress.report(get_spectrum_first_scan(spectrum), charge);
 
    // with the target database decide how many peptide mods to use
    MATCH_COLLECTION_T* target_psms = new_empty_match_collection(is_decoy); 
    int max_pep_mods = mpsm_search_pep_mods( target_psms, 
                                        is_decoy,   
                                        index,       
                                        database, 
                                        spectrum, 
                                        charge,
                                        peptide_mods, 
                                        num_peptide_mods,
                                        compute_pvalues); 
 

    // are there any matches?
    if( get_match_collection_match_total(target_psms) == 0 ){
      // don't print and don't search decoys
      carp(CARP_WARNING, "No matches found for spectrum %i, charge %i",
           get_spectrum_first_scan(spectrum), charge);
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

      mpsm_search_pep_mods(decoy_psms, 
                      is_decoy,   
                      index, 
                      database, 
                      spectrum, 
                      charge, 
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
  search_for_mpsms(spsm_map, mpsm_map);
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

  for (MPSM_ChargeMap::iterator iter = spsm_map.begin();
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
  free_spectrum_collection(spectra);
  free_parameters();

  exit(0);
}// end main
#else // SEARCH_ENABLED not defined
int mpsm_search_main(int argc, char **argv){
  (void) argc;
  (void) argv;
  fputs(
    "You are using the open source version of Crux. Due to intellectual\n"
    "property issues, we are unable to provide database search functionality\n"
    "in this version. To obtain a licence for the full functional version of\n"
    "Crux that includes the database search tools, please visit the following URL:\n"
    "\nhttp://depts.washington.edu/ventures/UW_Technology/Express_Licenses/crux.php\n",
    stderr
  );
  return 1;
}
#endif // SEARCH_ENABLED



/* Private function definitions */

/**
 * \brief Look at matches and search parameters to determine if a
 * sufficient number PSMs have been found.  Returns TRUE if the
 * maximum number of modifications per peptide have been considered.
 * In the future, implement and option and test for a minimum score.
 * \returns TRUE if no more PSMs need be searched.
 */
BOOLEAN_T mpsm_is_search_complete(MATCH_COLLECTION_T* matches, 
                             int mods_per_peptide){


  if( matches == NULL ){
    return FALSE;
  }

  // keep searching if no limits on how many mods per peptide
  if( get_int_parameter("max-mods") == MAX_PEPTIDE_LENGTH ){
    return FALSE;
  }
  // stop searching if at max mods per peptide
  if( mods_per_peptide == get_int_parameter("max-mods") ){ 
    return TRUE;
  }

  // test for minimun score found

  return FALSE;
  
}


/**
 * \brief Search the database OR index with up to num_peptide_mods from
 * the list for matches to the spectrum. 
 * Scored PSMs are added to the match_collection, possibly truncating
 * the collection and deleting existing matches in the collection.
 * After searching with each peptide mod, assess if there exists a
 * "good enough" match and end the search if there is, returning the
 * number of peptide mods that were searched.
 * \return The number of peptide mods searched.
 */
int mpsm_search_pep_mods(
  MATCH_COLLECTION_T* match_collection, ///< store PSMs here
  BOOLEAN_T is_decoy,   ///< generate decoy peptides from index/db
  INDEX_T* index,       ///< index to use for generating peptides
  DATABASE_T* database, ///< db to use for generating peptides
  SPECTRUM_T* spectrum, ///< spectrum to search
  int charge,           ///< seach spectrum at this charge state
  PEPTIDE_MOD_T** peptide_mods, ///< list of peptide mods to apply
  int num_peptide_mods, ///< how many p_mods to use from the list
  BOOLEAN_T store_scores///< save all scores for p-value estimation
){

  // set match_collection charge
  set_match_collection_charge(match_collection, charge);

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
      BOOLEAN_T passes = mpsm_is_search_complete(match_collection, cur_aa_mods);
      if( passes ){
        carp(CARP_DETAILED_DEBUG, 
             "Ending search with %i modifications per peptide", cur_aa_mods);
        break;
      }// else, search with more mods
      cur_aa_mods = this_aa_mods;
    }
      
    //TODO SJM:  Figure out why this code gives different results for the sequest 
    //smoke test (this was changed in Rev. 2006).
    //      20014c20014
    //< S     21134   21134   3       0.00    server  2140.03 0.00    0.00    213
    //---
    //> S     21134   21134   3       0.00    server  2140.03 0.00    0.00    207
    
    
    // get peptide iterator
    MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator =
      new_modified_peptides_iterator_from_mz(get_spectrum_precursor_mz(spectrum),
                                             charge,
                                             peptide_mod,
                                             is_decoy,
                                             index,
                                             database);
    
    
    // score peptides
  //FIX

 int num_matches_added = add_unscored_peptides(match_collection, spectrum, charge,
                                                peptide_iterator, is_decoy);

  // main scoring
  score_matches_one_spectrum(XCORR, match_collection, spectrum, charge, 
                             store_scores); 
  populate_match_rank_match_collection(match_collection, XCORR);

    free_modified_peptides_iterator(peptide_iterator);
    
  }//next peptide mod

  return mod_idx;
}

/**
 * Print the target and decoy match collections to their respective
 * target and decoy files.
 *
 * Three possibilities: 1. combine the target and all decoy
 * collections and print to target file.  2. print targets to target
 * file and combine all decoys and print to one decoy file.  3. print
 * each collection to a separate file.
 * Possible side effectos: Collections may be merged and re-ranked.
 */
/*
void print_spectrum_matches(
  OutputFiles& output_files,       
  MATCH_COLLECTION_T* target_psms, 
  MATCH_COLLECTION_T** decoy_psms,
  int num_decoy_collections,
  SPECTRUM_T* spectrum,             
  BOOLEAN_T combine_target_decoy,
  int num_decoy_files
                   ){

  // now print matches to one, two or several files
  if( combine_target_decoy == TRUE ){
    // merge all collections
    MATCH_COLLECTION_T* all_psms = target_psms;
    for(int decoy_idx = 0; decoy_idx < num_decoy_collections; decoy_idx++){
      merge_match_collections(decoy_psms[decoy_idx], all_psms);
    }
    
    // sort and rank
    if( get_match_collection_scored_type(all_psms, SP) == TRUE ){
      populate_match_rank_match_collection(all_psms, SP);
    }
    populate_match_rank_match_collection(all_psms, XCORR);
    
    output_files.writeMatches(all_psms, // target matches
                              NULL,     // decoy matches
                              0,        // num decoys
                              XCORR, spectrum); 
    
  }else{ // targets and decoys in separate files
    
    // if decoys in one file
    if( num_decoy_files == 1 ){
      // merge decoys
      MATCH_COLLECTION_T* merged_decoy_psms = decoy_psms[0];
      for(int decoy_idx = 1; decoy_idx < num_decoy_collections; decoy_idx++){
        merge_match_collections(decoy_psms[decoy_idx],
                                merged_decoy_psms);
      }
      
      // sort and rank
      if( get_match_collection_scored_type(merged_decoy_psms, SP) == TRUE ){
        populate_match_rank_match_collection(merged_decoy_psms, SP);
      }
      populate_match_rank_match_collection(merged_decoy_psms, XCORR);
      
      output_files.writeMatches(target_psms, &merged_decoy_psms, 
                                1, // num decoys
                                XCORR, spectrum);
      
    }else{
      // already sorted and ranked
      output_files.writeMatches(target_psms, decoy_psms, 
                                num_decoy_collections, XCORR, spectrum);
    }
  }
}
*/
// TODO this should be in match_collection
/**
 * Search the given database or index using shuffled peptides and the
 * spectrum/charge in the target psm match collection.  Add those
 * scores to the target psm match collection for use in weibull
 * parameter estimation but do not save the matches.  Repeat the
 * search with all peptide mods in the list.
 */
/*
void add_decoy_scores(
  MATCH_COLLECTION_T* target_psms, ///< add scores to these matches
  SPECTRUM_T* spectrum, ///<
  int charge, ///< 
  INDEX_T* index, ///< search this index if not null
  DATABASE_T* database, ///< search this database if not null
  PEPTIDE_MOD_T** peptide_mods, ///< list of peptide mods to search
  int num_peptide_mods ///< number of mods in the above array
){

  int mod_idx = 0;
  // for each peptide mod in the list
  for(mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++){
    MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator = 
      new_modified_peptides_iterator_from_mz(
                                          get_spectrum_precursor_mz(spectrum),
                                          charge,
                                          peptide_mods[mod_idx],
                                          index,
                                          database);
    add_decoy_scores_match_collection(target_psms, 
                                      spectrum, 
                                      charge, 
                                      peptide_iterator);  
  }


}
*/


bool passRTimeThreshold(MPSM_Match& match,
  FLOAT_T rtime_max_diff) {

  if (rtime_threshold) {
    //cerr <<"Testing rtime"<<endl;
    ChargeIndex& charge_index = match.getChargeIndex();

    double fdiff = fabs(rtime_max_diff);
    bool ans = false;
    //all +2
    if (charge_index.numCharge(2) == match.numMatches()) {
      //cerr <<"Done. All +2"<<endl;
      ans = fdiff <= rtime_all2_threshold; //8.4;
    } else if (charge_index.numCharge(3) == match.numMatches()) {
      //cerr <<"Done .All +3"<<endl;
      ans = fdiff <= rtime_all3_threshold; //14.1;
    } else {
      //cerr <<"Done. Mixture"<<endl;
      ans = fdiff <= rtime_default_threshold; //14.7;
    }
    return ans;
  } else {
    return true;
  }
}

RetentionPredictor* rtime_predictor = NULL;


BOOLEAN_T extendMatch(MPSM_Match& orig_mpsm, 
  MPSM_MatchCollection& spsm_matches,
  MPSM_ChargeMap& new_mpsms_matches,
  vector<set<string> >& visited,
  int match_collection_idx) {
  
  //cerr <<"extendMatch: start"<<endl;
  BOOLEAN_T match_added = FALSE;
  
  if (rtime_predictor == NULL) {
    rtime_predictor = RetentionPredictor::createRetentionPredictor();
  }
  set<string>& visited_here = visited.at(match_collection_idx);
  int n = spsm_matches.numMatches();

  if (get_boolean_parameter("mpsm-search-old-way")) {
    int top_n = get_int_parameter("mpsm-top-n");
    if (top_n != -1) {
      n = min(spsm_matches.numMatches(), top_n);
    }
  }

  for (int idx=0;idx < n; idx++) {
    MPSM_Match new_match(orig_mpsm);
    if (new_match.addMatch(spsm_matches.getMatch(idx).getMatch(0))) {

      if (visited_here.find(new_match.getString()) == visited_here.end()) {
        visited_here.insert(new_match.getString());
        FLOAT_T rtime_max_diff = rtime_predictor -> calcMaxDiff(new_match);

        if (passRTimeThreshold(new_match, rtime_max_diff)) {
          new_match.setRTimeMaxDiff(rtime_max_diff);
          new_mpsms_matches.insert(new_match, match_collection_idx);
          match_added = TRUE;
        }
      }
    }
  }
  //cerr <<"extendMatch:done. "<< match_added <<endl;
  return match_added;
}

BOOLEAN_T extendChargeMap(MPSM_ChargeMap& spsm_map,
                     MPSM_ChargeMap& current_mpsm_map,
                     int mpsm_level) {


  int top_n = get_int_parameter("mpsm-top-n");

  MPSM_ChargeMap::iterator map_iter;
  MPSM_ChargeMap::iterator map_iter2;

  set<ChargeIndex> charge_extended;

  vector<set<string> > visited;
  set<string> visited_target;
  visited.push_back(visited_target);
  for (int idx=0;idx < 3/*TODO - fix*/; idx++) {
    set<string> current_decoy;
    visited.push_back(current_decoy);
  }

  for (map_iter = current_mpsm_map.begin();
        map_iter != current_mpsm_map.end();
        ++map_iter) {

  

    ChargeIndex charge_index = map_iter -> first;

    carp(CARP_DEBUG,"mpsm_level:%d charge size:%d",mpsm_level, charge_index.size());
    if (charge_index.size() != mpsm_level) {
      continue;
    }
    
    if (charge_extended.find(charge_index) != charge_extended.end()) {
      continue;
    }

    //cout <<"Extending charge: "<<charge_index<<endl;
    charge_extended.insert(charge_index);



    vector<MPSM_MatchCollection>& match_collections = map_iter -> second;

 

    //start with targets.
    //Take the top candidate and consider matching it with every one else.
    MPSM_MatchCollection& target_collection = match_collections.at(0);
    
    if (get_boolean_parameter("mpsm-do-sort")) {
      target_collection.sortByScore(XCORR);
    }

    int current_top_n;
    if (top_n == -1) {
      current_top_n = target_collection.numMatches();
    } else {
      current_top_n = min(target_collection.numMatches(), top_n);
    }

    for (int current_index = 0;current_index < current_top_n; current_index++) {
      MPSM_Match& current_mpsm_target = target_collection.getMatch(current_index);
      carp(CARP_DEBUG,"searching %d of %d", current_index+1, current_top_n);

      //Loop through spsms by charge, and extend the current mpsm target using targets.
      for (map_iter2 = spsm_map.begin();
        map_iter2 != spsm_map.end();
        ++map_iter2) {

        ChargeIndex spsm_charge_index = map_iter2 -> first;

        ChargeIndex new_charge_index(charge_index);
        new_charge_index.add(spsm_charge_index);
        //cerr<<"Finding matches for:"<<new_charge_index<<endl;

        current_mpsm_map.insert(new_charge_index);

        vector<MPSM_MatchCollection>& spsm_match_collections = map_iter2 -> second;
        MPSM_MatchCollection& spsm_target_match_collection = spsm_match_collections.at(0);
      
        BOOLEAN_T matches_added = extendMatch(current_mpsm_target, 
          spsm_target_match_collection, 
          current_mpsm_map,
          visited,
          0);
        //cerr<<"Extending Target:"<<matches_added<<endl;
         //if there are some matches found, then search the decoys.
  
        if (!matches_added) {
          continue; //search the next charge state extension.
        }
      
        for (int idx=0;idx < match_collections.size()-1;idx++) {

          if (get_boolean_parameter("mpsm-do-sort")) {
            match_collections.at(idx+1).sortByScore(XCORR);
          }
        }
      
        for (int decoy_idx=1;decoy_idx < spsm_match_collections.size();decoy_idx++) {
          MPSM_MatchCollection& spsm_decoy_match_collection = spsm_match_collections.at(decoy_idx);
          if (get_boolean_parameter("mpsm-do-sort")) {
            spsm_decoy_match_collection.sortByScore(XCORR);
          }  
          if (get_boolean_parameter("mpsm-decoy-decoy")) {
            //decoy-decoy
            MPSM_MatchCollection& current_decoy_collection = match_collections.at(decoy_idx);
            MPSM_Match& current_mpsm_decoy = current_decoy_collection.getMatch(current_index);

            matches_added = extendMatch(current_mpsm_decoy,
              spsm_decoy_match_collection,
              current_mpsm_map,
              visited,
              decoy_idx);
            //cerr<<"Extending decoy("<<decoy_idx<<":"<<matches_added<<endl;

          } else {
            //target-decoy
            
            extendMatch(current_mpsm_target,
              spsm_decoy_match_collection,
              current_mpsm_map,
              visited,
              decoy_idx);
          }
       
        } /* for (decoy_idx)*/
      } /* map_iter2++ */
    } /* current_index++ */
  } /* map_iter++ */
  if (get_boolean_parameter("mpsm-do-sort")) {
    //cerr<<"Done extending (sorting)"<<endl;
    current_mpsm_map.sortMatches(XCORR);
  }
}


void search_for_mpsms(MPSM_ChargeMap& charge_spsm_map, 
  MPSM_ChargeMap& charge_mpsm_map) {

  charge_mpsm_map.clear();

  int max_peptides = get_int_parameter("mpsm-max-peptides");
 
  if (max_peptides < 2) return;
  
  charge_mpsm_map = charge_spsm_map;

  for (int npeptides = 2 ; npeptides <= max_peptides ; npeptides++) {
    //create new mpsm matches by extending the current
    //mpsm map by one using the spsm map
    MPSM_ChargeMap new_mpsm_map;
    BOOLEAN_T added = extendChargeMap(charge_spsm_map,
        charge_mpsm_map,
        npeptides-1);
    
    if (!added)
      break;
  }
  //cerr<<"Done searching for mpsms for this scan.."<<endl;
}

