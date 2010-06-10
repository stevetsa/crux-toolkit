/**
 * \file crux-main.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: November 24, 2008
 * \brief The starting point for the main crux program.
 *
 * Usage is "crux [command] [options] [arguments]" where command
 * is one of the primary crux commands.
 **/

#include "crux-main.h"

const char* usage_str = "Usage: crux <command> [options] <argument>\n"
"\n"
"Crux supports the following commands:\n"
"  create-index        Create an index for all peptides in a fasta file.\n"
"  search-for-matches  Search a collection of spectra against a sequence\n"
"                      database, returning a collection of peptide-spectrum\n"
"                      matches (PSMs) scored by XCorr.\n"
"  sequest-search      Similar to search-for-matches but use Sp as a \n"
"                      preliminary score followed by XCorr.\n"
"  compute-q-values    Assign a q-value, which is a statistical confidence\n"
"                      measure that accounts for multiple testing, to each\n"
"                      PSM in a given set.\n" 
"  percolator          Analyze a collection of PSMs to target and decoy\n"
"                      sequences using the percolator algorithm.\n"
"  q-ranker            Analyze a collection of PSMs using the Q-ranker\n"
"                      algorithm.\n"
"  print-processed-spectra\n"
"                      Write a new ms2 file with all of the same spectra\n"
"                      with only the peaks used for computing xcorr.\n"
"  search-for-xlinks   Search a collection of spectra against a sequence\n"
"                      database, returning a collection of matches\n"
"                      corresponding to linear and cross-linked peptides\n"
"                      scored by XCorr.\n"
"  version             Print the Crux version number to standard output,\n"
"                      then exit.\n"
"\n"
"Options and arguments are specific to each command. Type 'crux <command>'\n"
"for details.\n"
; 

/**
 * \brief Takes a directory containing PSM files and a protein index
 * and analyzes the PSMs using compute-q-values, percolator or q-ranker.
 */
static void analyze_matches_main(
  COMMAND_T command,
  int argc,
  char** argv
){

  // Define optional command line arguments.
  const char* qvalue_option_list[] = {
    "verbosity",
    "parameter-file",
    "overwrite",
    "output-dir",
    "fileroot"
  };
  int qvalue_num_options = 5;
  const char* percolator_option_list[] = {
    "verbosity",
    "parameter-file",
    "fileroot",
    "feature-file",
    "output-dir",
    "overwrite"
  };
  int percolator_num_options = 6;
  const char* qranker_option_list[] = {
    "verbosity",
    "parameter-file",
    "fileroot",
    "feature-file",
    "output-dir",
    "overwrite",
  };
  int qranker_num_options = 6;

  // Define required command line arguments.
  const char* argument_list[] = {
    "protein input",
  };
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  // Do some common initialization stuff.
  switch(command) {
  case QVALUE_COMMAND:
    initialize_run(command, argument_list, num_arguments,
		   qvalue_option_list, qvalue_num_options, argc, argv);
    break;
  case PERCOLATOR_COMMAND:
    initialize_run(command, argument_list, num_arguments,
		   percolator_option_list, percolator_num_options, argc, argv);
    break;
  case QRANKER_COMMAND:
    initialize_run(command, argument_list, num_arguments,
		   qranker_option_list, qranker_num_options, argc, argv);
    break;
  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;
  }

  // Get arguments
  char* psm_dir = get_string_parameter("output-dir");
  char* protein_input_name = get_string_parameter("protein input");

  // Prepare the output files.
  OutputFiles output(command);
  output.writeHeaders();

  // Perform the analysis.
  MATCH_COLLECTION_T* match_collection = NULL;
  switch(command) {
  case QVALUE_COMMAND:
    match_collection = run_qvalue(psm_dir,
				  protein_input_name);
    break;
  case PERCOLATOR_COMMAND:
    match_collection = run_percolator(psm_dir,
				      protein_input_name,
				      output);
    break;
  case QRANKER_COMMAND:
    match_collection = run_qranker(psm_dir,
				   protein_input_name,
				   output);
    break;
  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;
  }

  carp(CARP_INFO, "Outputting matches.");
  output.writeMatches(match_collection);

  // MEMLEAK below causes seg fault (or used to)
  // free_match_collection(match_collection);

  // clean up
  free(psm_dir);
  free(protein_input_name);

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  switch(command) {
  case QVALUE_COMMAND:
    carp(CARP_INFO, "Finished crux compute-q-values.");
    break;
  case PERCOLATOR_COMMAND:
    carp(CARP_INFO, "Finished crux percolator.");
    break;
  case QRANKER_COMMAND:
    carp(CARP_INFO, "Finished crux q-ranker.");
    break;
  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;
  }
}



int main(int argc, char** argv){

  // check the syntax for crux <operation>
  if( argc < 2 ){
    carp(CARP_FATAL, usage_str);
  }

  // determine the operation
  char* op_string = argv[1];
  COMMAND_T command = string_to_command_type(op_string);

  // call the appropriate function 
  // passing the command line minus the first token ('crux')
  switch(command){
  case INDEX_COMMAND:
    create_index_main(argc-1, argv+1);
    break;

  case SEARCH_COMMAND:
    search_main(argc-1, argv+1);
    break;

  case SEQUEST_COMMAND:
    sequest_search_main(argc-1, argv+1);
    break;

  case PROCESS_SPEC_COMMAND:
    print_processed_spectra_main(argc-1, argv+1);
    break;
  
  case XLINK_SEARCH_COMMAND:
    xlink_search_main(argc-1, argv+1);
    break;

  case QVALUE_COMMAND:
    analyze_matches_main(QVALUE_COMMAND, argc-1, argv+1);
    break;

  case QRANKER_COMMAND:
    analyze_matches_main(QRANKER_COMMAND, argc-1, argv+1);
    break;

  case PERCOLATOR_COMMAND:
    analyze_matches_main(PERCOLATOR_COMMAND, argc-1, argv+1);
    break;

  case VERSION_COMMAND:
    printf("Crux version %s\n", VERSION);
    break;    

  case INVALID_COMMAND:
    carp(CARP_FATAL, "Invalid command '%s'\n%s", op_string, usage_str);
    break;

  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;

  }

  exit (0);
}// end main















