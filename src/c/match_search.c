/**
 * \file match_search.c
 * BASED ON: original_match_search.c
 * DATE: Aug 19, 2008
 * AUTHOR: Barbara Frewen
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
#include "carp.h"
#include "parameter.h"
#include "spectrum_collection.h"
#include "match_collection.h"
#include "search_loop.h"

#define NUM_SEARCH_OPTIONS 16
#define NUM_SEARCH_ARGS 2
#define PARAM_ESTIMATION_SAMPLE_COUNT 500

/* Private functions */
int prepare_protein_input(char* input_file, 
                          INDEX_T** index, 
                          DATABASE_T** database);

// C interface to OutputFiles
void* construct_output_files(int num_decoys, int num_proteins);
void close_output_files(void* output_files_param, int num_spectra);

//int main(int argc, char** argv){
int search_main(int argc, char** argv){

  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);

  /* Define optional command line arguments */
  int num_options = NUM_SEARCH_OPTIONS;
  char* option_list[NUM_SEARCH_OPTIONS] = {
    "verbosity",
    "version",
    "parameter-file",
    "write-parameter-file",
    "overwrite",
    "use-index",
    /*
    "prelim-score-type",
    "score-type",
    */
    "compute-p-values",
    "spectrum-min-mass",
    "spectrum-max-mass",
    "spectrum-charge",
    "match-output-folder",
    "output-mode",
    "sqt-output-file",
    "tab-output-file",
    "decoy-sqt-output-file",
    "number-decoy-set"
  };

  /* Define required command line arguments */
  int num_arguments = NUM_SEARCH_ARGS;
  char* argument_list[NUM_SEARCH_ARGS] = {"ms2 file", "protein input"};

  /* Initialize parameter.c and set default values*/
  initialize_parameters();

  /* Define optional and required arguments */
  select_cmd_line_options(option_list, num_options);
  select_cmd_line_arguments(argument_list, num_arguments);

  /* Parse the command line, including optional params file
     Includes syntax, type, and bounds checking, dies on error */
  parse_cmd_line_into_params_hash(argc, argv, "crux search-for-matches");

  /* Set verbosity */
  //verbosity = get_int_parameter("verbosity");
  //set_verbosity_level(verbosity);

  /* Set seed for random number generation */
  if(strcmp(get_string_parameter_pointer("seed"), "time")== 0){
    time_t seconds; // use current time to seed
    time(&seconds); // Get value from sys clock and set seconds variable.
    srand((unsigned int) seconds); // Convert seconds to a unsigned int
  }
  else{
    srand((unsigned int)atoi(get_string_parameter_pointer("seed")));
  }
  
  carp(CARP_INFO, "Beginning crux search-for-matches");

  /* Get input: ms2 file */
  char* ms2_file = get_string_parameter_pointer("ms2 file");

  // open ms2 file
  SPECTRUM_COLLECTION_T* spectra = new_spectrum_collection(ms2_file);

  // parse the ms2 file for spectra
  carp(CARP_INFO, "Reading in ms2 file %s", ms2_file);
  if(!parse_spectrum_collection(spectra)){
    carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2_file);
    free_spectrum_collection(spectra);
    exit(1);
  }
  
  carp(CARP_DEBUG, "There were %i spectra found in the ms2 file",
       get_spectrum_collection_num_spectra(spectra));

  /* Get input: protein file */
  //char* input_file = get_string_parameter_pointer("protein input");
  char* input_file = get_string_parameter("protein input");

  /* Prepare input, fasta or index */
  INDEX_T* index = NULL;
  DATABASE_T* database = NULL;
  int num_proteins = prepare_protein_input(input_file, &index, &database); 
  free(input_file);

  carp(CARP_DEBUG, "Found %i proteins", num_proteins);
  if( num_proteins == 0 ){
    carp(CARP_FATAL, "No proteins were found in the protein source.");
    exit(1);
  }
  
  /* Prepare output files */

  // flags and counters for loop
  int num_decoys = get_int_parameter("number-decoy-set");
  void* output_files = construct_output_files(num_decoys, num_proteins);
  // eventually change to: OutputFiles output_files(num_decoys);

  /* Perform search: loop over spectra*/

  // create spectrum iterator
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* spectrum_iterator = 
    new_filtered_spectrum_charge_iterator(spectra);

  // get search parameters for match_collection
  BOOLEAN_T compute_pvalues = get_boolean_parameter("compute-p-values");
  int sample_count = (compute_pvalues) ? PARAM_ESTIMATION_SAMPLE_COUNT : 0;
  BOOLEAN_T combine_target_decoy = get_boolean_parameter("tdc");

  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );
  // for estimating params for p-values, randomly select a total of 
  //    sample_count matches, a constant fraction from each peptide mod
  int sample_per_pep_mod =  sample_count / num_peptide_mods;
  carp(CARP_DEBUG, "Got %d peptide mods, sample %i per", 
       num_peptide_mods, sample_per_pep_mod);

  int num_spectra = search_loop(spectrum_iterator,
				combine_target_decoy,
				num_peptide_mods,
				peptide_mods,
				database,
				index,
				sample_per_pep_mod,
				compute_pvalues,
				output_files,
				num_decoys);

  close_output_files(output_files, num_spectra);
  // Eventually: output_files.Close(num_spectra);

  // clean up memory

  carp(CARP_INFO, "Finished crux-search-for-matches");
  exit(0);
}// end main




/* Private function definitions */
/**
 * \brief Open either the index or fasta file and prepare it for
 * searching.  Die if the input file cannot be found or read.
 * \returns the number of proteins in the file/index
 */
int prepare_protein_input(char* input_file, 
                          INDEX_T** index, 
                          DATABASE_T** database){

  int num_proteins = 0;
  BOOLEAN_T use_index = get_boolean_parameter("use-index");

  if (use_index == TRUE){
    carp(CARP_INFO, "Preparing protein index %s", input_file);
    *index = new_index_from_disk(input_file);

    if (index == NULL){
      carp(CARP_FATAL, "Could not create index from disk for %s", input_file);
      exit(1);
    }
    num_proteins = get_index_num_proteins(*index);

  } else {
    carp(CARP_INFO, "Preparing protein fasta file %s", input_file);
    *database = new_database(input_file, FALSE);         
    if( database == NULL ){
      carp(CARP_FATAL, "Could not create protein database");
      exit(1);
    } 

    if(!parse_database(*database)){
      carp(CARP_FATAL, "Error with protein input");
      exit(1);
    } 
    num_proteins = get_database_num_proteins(*database);
  }
  return num_proteins;
}



