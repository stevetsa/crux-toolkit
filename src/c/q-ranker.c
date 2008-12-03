/**
 * \file q-ranker.c
 */
/*
 * AUTHOR: Marina Spivak
 * CREATE DATE: November 24, 2008
 * DESCRIPTION: ...
 * 
 * $Revision: 1.1.2.1 $
 ****************************************************************************/
#include "q-value.h"

#define NUM_QRANKER_OPTIONS 7
#define NUM_QRANKER_ARGUMENTS 2

/* Private function declarations. */
// if there are any

/**
 * \brief One of the commands for crux.  Takes in a directory
 * containing binary psm files and a protein source (index or fasta
 * file) and ... 
 */
int qranker_main(int argc, char** argv){

  /* Define command line arguments */
  int num_options = NUM_QRANKER_OPTIONS;
  char* option_list[NUM_QRANKER_OPTIONS] = {
    "version",                    // print version number
    "verbosity",                  // set level of output
    "parameter-file",             // use settings in file
    "write-parameter-file",       // write a file with settings
    "use-index",                  // index vs fasta file for proteins
    "overwrite",                  // clobber existing files?
    "sqt-output-file"             // name of output file
    // add more options here
  };

  int num_arguments = NUM_QRANKER_ARGUMENTS;
  char* argument_list[NUM_QRANKER_ARGUMENTS] = {
    "psm-folder",
    "protein input",
  };

  /* for debugging handling of parameters*/
  set_verbosity_level(CARP_ERROR);

  /* Set up parameters and set defaults in parameter.c */
  initialize_parameters();

  /* Define optional and required arguments in parameter.c */
  select_cmd_line_options(option_list, num_options );
  select_cmd_line_arguments(argument_list, num_arguments);

  /* Parse the command line and optional paramter file
     does sytnax, type, and bounds checking and dies on error */
  parse_cmd_line_into_params_hash(argc, argv, "crux q-ranker");

  /* Get arguments */
  char* psm_directory = get_string_parameter("psm-folder");
  char* protein_input_name = get_string_parameter("protein input");

  /* Perform the analysis */

  carp(CARP_INFO, "Running q-ranker");

  /*  The basic steps are to read all of the matches from the .csm
      files into some data structure that q-ranker can read, pass it
      to q-ranker, do magic, turn the results into a match_collection,
      and print the results. 
   */

  // this will iterate over all .csm files in the input folder
 MATCH_COLLECTION_ITERATOR_T* csm_file_iterator =
    new_match_collection_iterator(psm_directory, protein_input_name);

 while( match_collection_iterator_has_next(csm_file_iterator) ){

   // turn the file into a match_collection
   MATCH_COLLECTION_T* match_collection = 
     match_collection_iterator_next(csm_file_iterator);

   // this will iterate over all matches in the file/collection
   MATCH_ITERATOR_T* match_iterator = new_match_iterator(match_collection, 
                                                         XCORR,
                                                         FALSE);// don't sort

   while(match_iterator_has_next(match_iterator)){
     //MATCH_T* cur_match = match_iterator_next(match_iterator);

     // get information from match and put it in a data structure for q-ranker

   }// next match in file/collection

 }// next csm file

 // run q-ranker
 // turn results back into matches

  carp(CARP_INFO, "Outputting matches.");
  //print_sqt_file(match_collection, scorer_type, second_scorer_type);

  // clean up
  free(protein_input_name);

  carp(CARP_INFO, "crux q-ranker finished.");
  exit(0);
}

/*  ****************** Subroutines ****************/

/*
 */
void print_sqt_file_q(
  MATCH_COLLECTION_T* match_collection,
  SCORER_TYPE_T scorer,
  SCORER_TYPE_T second_scorer
  ){

  // get filename and open file
  char* sqt_filename = get_string_parameter("sqt-output-file");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  FILE* sqt_file = create_file_in_path( sqt_filename, NULL, overwrite );

  // print header
  int num_proteins = get_match_collection_num_proteins(match_collection);
  print_sqt_header( sqt_file, "target", num_proteins, TRUE);

  ALGORITHM_TYPE_T algorithm_type = get_algorithm_type_parameter("algorithm");
  char algorithm_str[64];
  algorithm_type_to_string(algorithm_type, algorithm_str);

  fprintf(sqt_file, "H\tComment\tmatches analyzed by %s\n", algorithm_str);

  // get match iterator sorted by spectrum
  MATCH_ITERATOR_T* match_iterator = 
    new_match_iterator_spectrum_sorted(match_collection, scorer);

  // print each spectrum only once, keep track of which was last printed
  int cur_spectrum_num = -1;
  int cur_charge = 0;
  int match_counter = 0;
  int max_matches = get_int_parameter("max-sqt-result");

  // for all matches
  while( match_iterator_has_next(match_iterator) ){

    // get match and spectrum
    MATCH_T* match = match_iterator_next(match_iterator);
    SPECTRUM_T* spectrum = get_match_spectrum(match);
    int this_spectrum_num = get_spectrum_first_scan(spectrum);
    int charge = get_match_charge(match);

    carp(CARP_DETAILED_DEBUG, 
         "SQT printing scan %i (current %i), charge %i (current %i)", 
         this_spectrum_num, cur_spectrum_num, charge, cur_charge);

    // if this spectrum has not been printed...
    if( cur_spectrum_num != this_spectrum_num
        || cur_charge != charge){

      carp(CARP_DETAILED_DEBUG, "Printing new S line");
      // print S line to sqt file
      cur_spectrum_num = this_spectrum_num;
      cur_charge = charge;
      int num_peptides = get_match_ln_experiment_size(match);
      num_peptides = expf(num_peptides);

      print_spectrum_sqt(spectrum, sqt_file, num_peptides, charge);

      // print match to sqt file
      print_match_sqt(match, sqt_file, scorer, second_scorer);
      match_counter = 1;
    }
    // if this spectrum has been printed
    else{  
      if( match_counter < max_matches ){
        print_match_sqt(match, sqt_file, scorer, second_scorer);
        match_counter++;
      }
    }

  }// next match
  free_match_iterator(match_iterator);
  free(sqt_filename);

}






/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
