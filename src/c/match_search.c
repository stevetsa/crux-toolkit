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
void open_output_files(FILE*** binary_filehandle_array, 
                       FILE** sqt_filehandle,
                       FILE** decoy_sqt_filehandle,
                       FILE** tab_file,
                       FILE** decoy_tab_file);


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

  FILE** psm_file_array = NULL; //file handle array
  FILE* sqt_file = NULL;
  FILE* decoy_sqt_file  = NULL;
  FILE* tab_file = NULL;
  FILE* decoy_tab_file  = NULL;

  open_output_files(
    &psm_file_array, 
    &sqt_file, 
    &decoy_sqt_file, 
    &tab_file,
    &decoy_tab_file
  );

  //print headers
  serialize_headers(psm_file_array);
  print_sqt_header(sqt_file, "target", num_proteins, FALSE);// !analyze-matches
  print_sqt_header(decoy_sqt_file, "decoy", num_proteins, FALSE);
  print_tab_header(tab_file);
  print_tab_header(decoy_tab_file);
  /* Perform search: loop over spectra*/

  // create spectrum iterator
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* spectrum_iterator = 
    new_filtered_spectrum_charge_iterator(spectra);

  // get search parameters for match_collection
  BOOLEAN_T compute_pvalues = get_boolean_parameter("compute-p-values");
  int sample_count = (compute_pvalues) ? PARAM_ESTIMATION_SAMPLE_COUNT : 0;
  BOOLEAN_T combine_target_decoy = get_boolean_parameter("tdc");

  // flags and counters for loop
  int num_decoys = get_int_parameter("number-decoy-set");

  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );
  // for estimating params for p-values, randomly select a total of 
  //    sample_count matches, a constant fraction from each peptide mod
  int sample_per_pep_mod =  sample_count / num_peptide_mods;
  carp(CARP_DEBUG, "Got %d peptide mods, sample %i per", 
       num_peptide_mods, sample_per_pep_mod);

  int spectrum_searches_counter = search_loop(spectrum_iterator,
					      combine_target_decoy,
					      num_peptide_mods,
					      peptide_mods,
					      database,
					      index,
					      sample_per_pep_mod,
					      compute_pvalues,
					      psm_file_array,
					      sqt_file,
					      decoy_sqt_file,
					      tab_file,
					      decoy_tab_file,
					      num_decoys);

  // fix headers in csm files
  int file_idx;
  for(file_idx=0; file_idx < num_decoys + 1; file_idx++){
    carp(CARP_DEBUG, "Changing csm header to have %i spectrum searches",
         spectrum_searches_counter);
    serialize_total_number_of_spectra(spectrum_searches_counter,
                                      psm_file_array[file_idx]);
  }
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

/**
 * \brief A private function for crux-search-for-matches to prepare
 * binary psm and text sqt files.
 *
 * Reads the --overwrite and --output-mode values from
 * parameter.c. Opens psm file(s) if requested, setting a given
 * pointer to the array of filehandles.  Opens sqt file(s) if
 * requested, setting the given pointers to each file handle.  If
 * binary files not requested, creates an array of NULL pointers.  If
 * sqt files not requested, sets given pointers to NULL. 
 *
 * \returns void.  Sets given arguments to newly created filehandles.
 */
void open_output_files(
  FILE*** psm_file_array, ///< put binary psm filehandles here -out
  FILE** sqt_file,        ///< put text sqt filehandle here -out
  FILE** decoy_sqt_file,  ///< put decoy sqt filehandle here -out
  FILE** tab_file,        ///< put text sqt filehandle here -out
  FILE** decoy_tab_file)  ///< put decoy sqt filehandle here -out
{
  char* match_output_folder = get_string_parameter("match-output-folder");
  MATCH_SEARCH_OUTPUT_MODE_T output_type = get_output_type_parameter(
                                                    "output-mode");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  carp(CARP_DEBUG, "The output type is %d (binary, sqt, tab, all)" \
       " and overwrite is '%d'", (int)output_type, (int)overwrite);


  // create binary psm files (allocate memory, even if not used)
  *psm_file_array = create_psm_files();

  if(output_type == SQT_OUTPUT || output_type == ALL_OUTPUT){

    //create sqt file handles
    carp(CARP_DEBUG, "Opening sqt files");
    char* sqt_filename = get_string_parameter_pointer("sqt-output-file");
    *sqt_file = create_file_in_path(sqt_filename, 
                                    match_output_folder, 
                                    overwrite);
    char* decoy_sqt_filename = get_string_parameter_pointer(
                                                    "decoy-sqt-output-file");
    if( get_int_parameter("number-decoy-set") > 0 ){
      *decoy_sqt_file = create_file_in_path(decoy_sqt_filename,
                                            match_output_folder,
                                            overwrite);
    }

    if(sqt_file == NULL || decoy_sqt_file == NULL){
      carp(CARP_DEBUG, "sqt file or decoy is null");
    }
  }

  if(output_type == TAB_OUTPUT || output_type == ALL_OUTPUT){

    //create sqt file handles
    carp(CARP_DEBUG, "Opening tab delimited files");
    char* tab_filename = get_string_parameter_pointer("tab-output-file");
    *tab_file = create_file_in_path(tab_filename, 
                                    match_output_folder, 
                                    overwrite);
    char* decoy_tab_filename = get_string_parameter_pointer(
                                                    "decoy-tab-output-file");
    if( get_int_parameter("number-decoy-set") > 0 ){
      *decoy_tab_file = create_file_in_path(decoy_tab_filename,
                                            match_output_folder,
                                            overwrite);
    }

    if(tab_file == NULL || decoy_tab_file == NULL){
      carp(CARP_DEBUG, "tab file or decoy tab file is null");
    }

  }

  free(match_output_folder);
  carp(CARP_DEBUG, "Finished opening output files");
}


