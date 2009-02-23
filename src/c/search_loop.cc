// search_loop moved out of match_search.c

extern "C" {
#include "carp.h"
#include "spectrum_collection.h"
#include "match_collection.h"
#include "search_loop.h"
}
#include "output_files.h"

BOOLEAN_T is_search_complete(MATCH_COLLECTION_T* matches, 
                             int mods_per_peptide);

extern "C" int search_loop(FILTERED_SPECTRUM_CHARGE_ITERATOR_T* spectrum_iterator,
			   BOOLEAN_T combine_target_decoy,
			   int num_peptide_mods,
			   PEPTIDE_MOD_T** peptide_mods,
			   DATABASE_T* database,
			   INDEX_T* index,
			   int sample_per_pep_mod,
			   BOOLEAN_T compute_pvalues,
			   void* output_files_param,
			   int num_decoys) {
  OutputFiles* output_files = (OutputFiles*) output_files_param;

  int mod_idx = 0;
  int spectrum_searches_counter = 0; //for psm file header, spec*charges

  // for each spectrum
  while(filtered_spectrum_charge_iterator_has_next(spectrum_iterator)){
    int charge = 0;
    SPECTRUM_T* spectrum = 
      filtered_spectrum_charge_iterator_next(spectrum_iterator, &charge);
    double mass = get_spectrum_neutral_mass(spectrum, charge);

    carp(CARP_DETAILED_INFO, 
         "Searching spectrum number %i, charge %i, search number %i",
         get_spectrum_first_scan(spectrum), charge,
         spectrum_searches_counter+1 );

    // with just the target database decide how many peptide mods to use
    // create an empty match collection
    MATCH_COLLECTION_T* match_collection = 
      new_empty_match_collection( FALSE ); // is decoy = false

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
        BOOLEAN_T passes = is_search_complete(match_collection, cur_aa_mods);
        if( passes ){
          carp(CARP_DETAILED_DEBUG, 
               "Ending search with %i modifications per peptide", cur_aa_mods);
          break;
        }// else, search with more mods
        cur_aa_mods = this_aa_mods;
      }

      // get peptide iterator
      MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator = 
        new_modified_peptides_iterator_from_mass(mass,
                                                 peptide_mod,
                                                 index,
                                                 database);
      // score peptides
      int added = add_matches(match_collection, 
                              spectrum, 
                              charge, 
                              peptide_iterator,
                              sample_per_pep_mod,
                              FALSE // is decoy
                              );

      carp(CARP_DEBUG, "Added %i matches", added);

      free_modified_peptides_iterator(peptide_iterator);

    }//next peptide mod

    // in case we searched all mods, do we need to assess again?

    // are there any matches?
    if( get_match_collection_match_total(match_collection) == 0 ){
      // don't print and don't search decoys
      carp(CARP_WARNING, "No matches found for spectrum %i, charge %i",
           get_spectrum_first_scan(spectrum), charge);
      free_match_collection(match_collection);
      continue; // next spectrum
    }
    
    // calculate p-values
    if( compute_pvalues == TRUE ){
      carp(CARP_DEBUG, "Estimating Weibull parameters.");
      if( estimate_weibull_parameters_from_sample_matches(match_collection,
                                                          spectrum,
                                                          charge) ){
        carp(CARP_DEBUG, "Calculating p-values.");
        compute_p_values(match_collection);
      }else{
        set_p_values_as_unscored(match_collection);
      }
    }

    if( combine_target_decoy == FALSE ){
      // print matches
      carp(CARP_DEBUG, "About to print target matches");
      output_files->PrintMatches(match_collection, spectrum);
      //does this free all the matches, all the spectra and all the peptides?
      free_match_collection(match_collection);
      match_collection = NULL;
    }
    // now score same number of mods for decoys
    int max_mods = mod_idx;

    // for num_decoys  and num_decoy_repeats
    int decoy_idx = 0;
    int repeat_idx = 0;
    int num_decoy_repeats = get_int_parameter("num-decoys-per-target");
    // decoy files is 0 but we want to do at least one decoy searc for tdc
    num_decoys = (combine_target_decoy) ? 1 : num_decoys;

    for(decoy_idx = 0; decoy_idx < num_decoys; decoy_idx++ ){
      carp(CARP_DETAILED_DEBUG, "Searching decoy %i", decoy_idx+1);

      if( combine_target_decoy == FALSE ){// create an empty match collection 
        match_collection = new_empty_match_collection( TRUE ); // is decoy
      }// else add to match_collection from above

      for(mod_idx = 0; mod_idx < max_mods; mod_idx++){
          
        // get peptide mod
        PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];
        
        // for multiple decoy searches written to one file, repeat here
        for(repeat_idx=0; repeat_idx < num_decoy_repeats; repeat_idx++){
          // get peptide iterator
          MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator = 
            new_modified_peptides_iterator_from_mass(mass,
                                                     peptide_mod,
                                                     index,
                                                     database);
          // score peptides
          int added = add_matches(match_collection, 
                                  spectrum, 
                                  charge, 
                                  peptide_iterator,
                                  0, // no sampling for param estimation
                                  TRUE);// is decoy
          carp(CARP_DEBUG, "Added %i matches", added);
          
          free_modified_peptides_iterator(peptide_iterator);
        }// next repeat
        
      }// last mod

      // print matches
      if( combine_target_decoy == FALSE ){ // print to decoy file
        carp(CARP_DEBUG, "About to print decoy matches");
        output_files->PrintMatches(match_collection, 
				   spectrum, 
				   true, // is decoy
				   1+decoy_idx,
				   decoy_idx == 0 // send_to_sqt_tab?
				   // (we only print first decoy to sqt)
				   );
        
        free_match_collection(match_collection);
      }else{ // print all to target file
	carp(CARP_DEBUG, "About to print target and decoy matches");
	output_files->PrintMatches(match_collection, spectrum);
      }

    }// last decoy

    spectrum_searches_counter++;

    // clean up
    //    free_match_collection(match_collection);
  }// next spectrum

  return spectrum_searches_counter;
}



/**
 * \brief Look at matches and search parameters to determine if a
 * sufficient number PSMs have been found.  Returns TRUE if the
 * maximum number of modifications per peptide have been considered.
 * In the future, implement and option and test for a minimum score.
 * \returns TRUE if no more PSMs need be searched.
 */
BOOLEAN_T is_search_complete(MATCH_COLLECTION_T* matches, 
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
