/*************************************************************************//**
 * \file q-value.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: Jan 03 2007
 * \brief  Given as input a directory containing binary psm files,
 * a protein database, and an optional parameter file, analyze the
 * matches (with percolator or q-value) and return scores indicating
 * how good the matches are. 
 *
 * Handles at most 4 files (target and decoy).  Expects psm files to
 * start with <fileroot>.se and 
 * end with the extension '.txt' and decoys to end with
 * '-decoy#.txt'.  Multiple target files in the given directory are
 * concatinated together and presumed to be non-overlaping parts of
 * the same ms2 file. 
 ****************************************************************************/
#include "q-value.h"
#include <map>

using namespace std;

static const int MAX_PSMS = 10000000;
// 14th decimal place
static const double EPSILON = 0.00000000000001;

/**
* Find the best-scoring match, by XCorr, for each peptide in a given
* collection.  Only consider the top-ranked PSM per spectrum.
*
* Results are stored in the given match collection.
*/
static void identify_best_psm_per_peptide
 (MATCH_COLLECTION_T* all_matches)
{
  /* Instantiate a hash table.  key = peptide; value = maximal xcorr
     for that peptide. */
  map<string, FLOAT_T> best_xcorr_per_peptide;

  // Store in the hash the best xcorr per peptide.
  MATCH_ITERATOR_T* match_iterator 
    = new_match_iterator(all_matches, XCORR, FALSE);
  while(match_iterator_has_next(match_iterator)){
    MATCH_T* match = match_iterator_next(match_iterator);

    // Skip matches that are not top-ranked.
    if (get_match_rank(match, XCORR) == 1) {
      char *peptide = get_match_sequence(match);
      FLOAT_T xcorr = get_match_score(match, XCORR);

      map<string, FLOAT_T>::iterator map_position 
	= best_xcorr_per_peptide.find(peptide);

      if (map_position == best_xcorr_per_peptide.end()) {
	best_xcorr_per_peptide[peptide] = xcorr;
      } else {
        if (map_position->second < xcorr) {
          best_xcorr_per_peptide[peptide] = xcorr;
        }
      }
      free(peptide);
    }
  }
  free_match_iterator(match_iterator);


  // Set the best_per_peptide Boolean in the match, based on the hash.
  match_iterator = new_match_iterator(all_matches, XCORR, FALSE);
  while(match_iterator_has_next(match_iterator)){
    MATCH_T* match = match_iterator_next(match_iterator);

     // Skip matches that are not top-ranked.
    if (get_match_rank(match, XCORR) == 1) {
      char* peptide = get_match_sequence(match);
      FLOAT_T xcorr = get_match_score(match, XCORR);

      map<string, FLOAT_T>::iterator map_position 
	= best_xcorr_per_peptide.find(peptide);

      if (map_position->second == xcorr) {
	set_best_per_peptide(match);

        // Prevent ties from causing two peptides to be best.
        best_xcorr_per_peptide[peptide] = HUGE_VAL;
     }

     free(peptide);
    }
  }

  free_match_iterator(match_iterator);
}

#ifdef OLD_VERSION
/**
 * Find the best-scoring match, by XCorr, for each peptide in a given
 * collection.  Only consider the top-ranked PSM per spectrum.
 *
 * Results are stored in the given match collection.
 */
static void identify_best_psm_per_peptide
  (MATCH_COLLECTION_T* all_matches)
{
  /* Instantiate a hash table.  key = peptide; value = maximal xcorr
     for that peptide. */
  HASH_T* best_xcorr_per_peptide 
    = new_hash(get_match_collection_match_total(all_matches));

  MATCH_ITERATOR_T* match_iterator 
    = new_match_iterator(all_matches, XCORR, FALSE);
  while(match_iterator_has_next(match_iterator)){
    MATCH_T* match = match_iterator_next(match_iterator);

    // Skip matches that are not top-ranked.
    if (get_match_rank(match, XCORR) == 1) {
      char* peptide = get_match_sequence(match);
      FLOAT_T xcorr = get_match_score(match, XCORR);

      FLOAT_T* best_xcorr = 
	(FLOAT_T*)get_hash_value(best_xcorr_per_peptide, peptide);
      if (best_xcorr == NULL) {
	add_hash(best_xcorr_per_peptide, peptide, (void*)&xcorr);
      } else if (*best_xcorr < xcorr) {
	update_hash_value(best_xcorr_per_peptide, peptide, (void*)&xcorr);
      }
    }
  }
  free_match_iterator(match_iterator);

  // Set the best_per_peptide Boolean in the match, based on the hash.
  match_iterator = new_match_iterator(all_matches, XCORR, FALSE);
  while(match_iterator_has_next(match_iterator)){
    MATCH_T* match = match_iterator_next(match_iterator);

    // Skip matches that are not top-ranked.
    if (get_match_rank(match, XCORR) == 1) {
      char* peptide = get_match_sequence(match);
      FLOAT_T xcorr = get_match_score(match, XCORR);
      FLOAT_T* stored_best = 
	(FLOAT_T*)get_hash_value(best_xcorr_per_peptide, peptide);

      if (*stored_best == xcorr) {
	set_best_per_peptide(match);

	// Prevent ties from causing two peptides to be best.
	FLOAT_T huge = HUGE_VAL;
	update_hash_value(best_xcorr_per_peptide, peptide, (void*)&huge);
      }
    }
    free_match_iterator(match_iterator);
  }
  free_hash(best_xcorr_per_peptide);
}
#endif

/**
 * Compare doubles
 */
static int compare_doubles_descending(
    const void *a,
    const void *b
    ){
  double temp = *((double *)a) - *((double *)b);
  if (temp > 0){
    return -1;
  } else if (temp < 0){
    return 1;
  } else {
    return 0;
  }
}

/**
 * The q-value is defined as the minimum FDR at which a given score is
 * deemed significant.  This function takes a list of FDRs in sorted
 * order and converts them into q-values.
 *
 * All calculations are done as negative logs.
 */
static void convert_fdr_to_qvalue 
  (FLOAT_T* qvalues,     ///< Come in as FDRs, go out as q-values.
   int      num_values)
{
  double max_log_qvalue = - BILLION;
  int idx;
  for (idx=num_values-1; idx >= 0; idx--){
    if (qvalues[idx] > max_log_qvalue){
      max_log_qvalue = qvalues[idx];
    } else { // current q-value is <= max q-value
      // set to max q-value so far
      qvalues[idx] = max_log_qvalue; 
    }
    carp(CARP_DETAILED_DEBUG, "qvalue[%i] = %.10f", idx, qvalues[idx]);
  }
}

/**
 * Store two parallel arrays of floats in a hash table.
 *
 * The new hash table must be freed by the caller.
 */
static HASH_T* store_arrays_as_hash
  (FLOAT_T* keys, 
   FLOAT_T* values,
   int      num_values
){

  HASH_T* return_value = new_hash(num_values);
  int idx;
  for (idx=0; idx < num_values; idx++){
    // QUESTION: How long should this string be?  And who will free it?
    char* key_string = (char*)mymalloc(30 * sizeof(char));
    // QUESTION: How much precision should I use to store this value?
    sprintf(key_string, "%g", keys[idx]);
    FLOAT_T value = values[idx];
    add_hash(return_value, key_string, (void*)&value);
  }
  return(return_value);
}

/**
 * Use the Benjamini-Hochberg procedure to convert a given set of
 * p-values into q-values.  
 *
 * Assumes that the input is an array of negative log p-values.  The
 * output array also contains negative log q-values.
 *
 * This function uses the command line parameter "pi-zero".
 */
static FLOAT_T* compute_qvalues_from_pvalues(
  FLOAT_T* pvalues, 
  int      num_pvals
){

  // sort the - log p-values in descending order
  qsort(pvalues, num_pvals, sizeof(double), compare_doubles_descending);

  // work in negative log space, since that is where p- and qvalues end up
  double log_num_psms = - log(num_pvals);
  double log_pi_0 = - log(get_double_parameter("pi-zero"));

  // convert the p-values into FDRs using Benjamini-Hochberg
  FLOAT_T* qvalues = (FLOAT_T*)mycalloc(num_pvals, sizeof(FLOAT_T));
  int idx;
  for (idx=0; idx < num_pvals; idx++){
    carp(CARP_DETAILED_DEBUG, "pvalue[%i] = %.10f", idx, pvalues[idx]);
    int pvalue_idx = idx + 1; // start counting pvalues at 1
    double log_pvalue = pvalues[idx];

    double log_qvalue = 
      log_pvalue + log_num_psms - (-log(pvalue_idx)) + log_pi_0;
    qvalues[idx] = log_qvalue;
    carp(CARP_DETAILED_DEBUG, "FDR[%i] = %.10f", idx, qvalues[idx]);
  }

  // convert the FDRs into q-values
  convert_fdr_to_qvalue(qvalues, num_pvals);

  return(qvalues);
}

/**
 * Given a hash table that maps from a score to its q-value, assign
 * q-values to all of the matches in a given collection.
 */
static void assign_qvalues(
  HASH_T* score_to_qvalue_hash,
  SCORER_TYPE_T score_type,
  MATCH_COLLECTION_T* all_matches
){

  // Iterate over the matches filling in the q-values
  MATCH_ITERATOR_T* match_iterator = new_match_iterator(all_matches, 
                                                        XCORR, FALSE);
  while(match_iterator_has_next(match_iterator)){
    MATCH_T* match = match_iterator_next(match_iterator);
    FLOAT_T score = get_match_score(match, score_type);

    // Retrieve the corresponding q-value.
    char score_string[30];
    sprintf(score_string, "%g", score);
    FLOAT_T* qvalue = (FLOAT_T*)get_hash_value(score_to_qvalue_hash, 
					       score_string);

    /* If we're given a base score, then store the q-value.  If we're
       given a q-value, then store the peptide-level q-value. */
    switch (score_type) {

    case XCORR:
      set_match_score(match, DECOY_XCORR_QVALUE, *qvalue);
      break;

    case DECOY_XCORR_QVALUE:
      set_match_score(match, DECOY_XCORR_PEPTIDE_QVALUE, *qvalue);
      break;

    case LOGP_BONF_WEIBULL_XCORR: 
      set_match_score(match, LOGP_QVALUE_WEIBULL_XCORR, *qvalue);
      break;

    case LOGP_QVALUE_WEIBULL_XCORR:
      set_match_score(match, LOGP_PEPTIDE_QVALUE_WEIBULL, *qvalue);
      break;

    case PERCOLATOR_SCORE:
      set_match_score(match, PERCOLATOR_QVALUE, *qvalue);
      break;

    case PERCOLATOR_QVALUE:
      set_match_score(match, PERCOLATOR_PEPTIDE_QVALUE, *qvalue);
      break;

    case QRANKER_SCORE:
      set_match_score(match, QRANKER_QVALUE, *qvalue);
      break;

    case QRANKER_QVALUE:
      set_match_score(match, QRANKER_PEPTIDE_QVALUE, *qvalue);
      break;

    // Should never reach this point.
    case SP: 
    case LOGP_WEIBULL_XCORR: 
    case DECOY_XCORR_PEPTIDE_QVALUE:
    case LOGP_PEPTIDE_QVALUE_WEIBULL:
    case PERCOLATOR_PEPTIDE_QVALUE:
    case QRANKER_PEPTIDE_QVALUE:
    case NUMBER_SCORER_TYPES:
      carp(CARP_FATAL, "Something is terribly wrong!");
    }
  }
}

/**
 * \brief Compute q-values from a given set of scores, using a second
 * set of scores as an empirical null.  Sorts the incoming target
 * scores and returns a corresponding list of q-values.
 */
static FLOAT_T* compute_decoy_qvalues(
  FLOAT_T* target_scores,
  int      num_targets,
  FLOAT_T* decoy_scores,
  int      num_decoys
){
  if ((num_targets == 0) || (num_decoys == 0)) {
    carp(CARP_FATAL, "Cannot compute q-values (%d targets, %d nulls).",
	 num_targets, num_decoys);
  }
  carp(CARP_DEBUG, "Computing decoy q-values.");

  // Sort both sets of scores.
  qsort(target_scores, num_targets, sizeof(FLOAT_T), compare_doubles_descending);
  qsort(decoy_scores, num_targets, sizeof(FLOAT_T), compare_doubles_descending);

  // Compute false discovery rate for each target score.
  FLOAT_T* qvalues = (FLOAT_T*)mycalloc(num_targets, sizeof(FLOAT_T));
  int decoy_idx = 0;
  int target_idx;
  for (target_idx = 0; target_idx < num_targets; target_idx++) {
    FLOAT_T target_score = target_scores[target_idx];

    // Find the index of the first decoy score 
    while ((decoy_idx < num_decoys) &&
	   (decoy_scores[decoy_idx] > target_score)) {
      decoy_idx++;
    }

    // FDR = #decoys / #targets
    FLOAT_T fdr = decoy_idx / target_idx;
    if( target_idx == 0 ){ fdr = 1.0; }

    qvalues[target_idx] = fdr;
  }

  // convert the FDRs into q-values
  convert_fdr_to_qvalue(qvalues, num_targets);

  return(qvalues);
}

/**
 * \brief Compute a q-values based on what is in the PSM files in the
 * directory.  Store q-values in the match collection returned.
 *
 * If p-values were computed, then perform Benjamini-Hochberg q-value
 * calculations. Otherwise, if decoys are present, then rank on xcorr
 * and compute empirical q-values based on the number of decoys and
 * targets above the score threshold.
 *
 * \returns a collection of target PSMs with one q-value in each
 * match.
 */
MATCH_COLLECTION_T* run_qvalue(
  char* input_directory, 
  char* fasta_file 
  ){

  int num_decoys = 0; // to be set by match_collection_iterator
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator =
    new_match_collection_iterator(input_directory, fasta_file, &num_decoys);
  if( num_decoys > 1 ){
    carp(CARP_FATAL, "Only one decoy file per target can be processed "
         "but %d were found.  Please move extra decoy files.", num_decoys);
  }

  // Create two match collections, for targets and decoys.
  MATCH_COLLECTION_T* target_matches = new_empty_match_collection(FALSE);
  MATCH_COLLECTION_T* decoy_matches = new_empty_match_collection(TRUE);
  set_match_collection_scored_type(target_matches, XCORR, TRUE);
  set_match_collection_scored_type(decoy_matches, XCORR, TRUE);

  // Did we find something from which to get q-values?
  BOOLEAN_T have_pvalues = FALSE;
  BOOLEAN_T have_decoys = FALSE;

  // Iterate over all match collections in this directory.
  while(match_collection_iterator_has_next(match_collection_iterator)){
    MATCH_COLLECTION_T* match_collection = 
      match_collection_iterator_next(match_collection_iterator);

    // Keep track of whether we got p-values.
    // N.B. Assumes that if one collection has p-values, they all do.
    have_pvalues = get_match_collection_scored_type(match_collection,
						    LOGP_BONF_WEIBULL_XCORR); 

    // Iterate, gathering matches into one or two collections.
    MATCH_ITERATOR_T* match_iterator =
      new_match_iterator(match_collection, XCORR, FALSE);
    while(match_iterator_has_next(match_iterator)){
      MATCH_T* match = match_iterator_next(match_iterator);

      // Only use top-ranked matches.
      if( get_match_rank(match, XCORR) != 1 ){
        continue;
      }

      if (get_match_null_peptide(match) == TRUE) {
	add_match_to_match_collection(decoy_matches, match);
	have_decoys = TRUE;
      } else {
	add_match_to_match_collection(target_matches, match);
      }
    }
    free_match_iterator(match_iterator);
  }
  free_match_collection_iterator(match_collection_iterator);

  // Compute q-values from p-values.
  FLOAT_T* pvalues = NULL; // N.B. Misnamed for decoy calculation.
  int num_pvals = get_match_collection_match_total(target_matches);
  FLOAT_T* qvalues = NULL;
  if (have_pvalues == TRUE) {
    carp(CARP_DEBUG, "There are %d PSMs for q-value computation.", num_pvals);
    set_match_collection_scored_type(target_matches, 
				     LOGP_BONF_WEIBULL_XCORR, 
				     TRUE);
    pvalues = extract_scores_match_collection(LOGP_BONF_WEIBULL_XCORR,
					      target_matches);
    qvalues = compute_qvalues_from_pvalues(pvalues, num_pvals);
  }

  // Compute q-values from the XCorr decoy distribution.
  else if (have_decoys == TRUE) {
    int num_decoys = get_match_collection_match_total(decoy_matches);
    carp(CARP_DEBUG,
	 "There are %d target and %d decoy PSMs for q-value computation.",
	 num_pvals, num_decoys);
    pvalues = extract_scores_match_collection(XCORR, target_matches);
    FLOAT_T* decoy_xcorrs 
      = extract_scores_match_collection(XCORR, decoy_matches);
    qvalues = compute_decoy_qvalues(pvalues, num_pvals, 
				    decoy_xcorrs, num_decoys);
    free(decoy_xcorrs);
  }

  // Fatal: Cannot compute q-values.
  else {
    carp(CARP_FATAL, "Cannot compute q-values without decoy PSMs or p-values.");
  }

  // Store p-values to q-values as a hash, and then assign them.
  HASH_T* qvalue_hash 
    = store_arrays_as_hash(pvalues, qvalues, num_pvals);
  assign_qvalues(qvalue_hash, XCORR, target_matches);
  free(pvalues);
  free(qvalues);
  free_hash(qvalue_hash);

  // Identify PSMs that are top-scoring per peptide.
  identify_best_psm_per_peptide(target_matches);

  // Compute peptide-level q-values.
  //  compute_decoy_q_values(all_matches/
  //			 TRUE); // Do peptide-level scoring.

  free_match_collection(decoy_matches);
  free_match_collection(target_matches);
  return(target_matches);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
