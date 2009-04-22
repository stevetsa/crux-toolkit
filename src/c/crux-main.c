/**
 * \file crux-main.c
 */
/*
 AUTHOR: Barbara Frewen
 CREATE DATE: November 24, 2008
 DESCRIPTION: The starting point for what were previously three
         separate protgrams--create-index, search-for-matches,
         analyze-matches.  Usage is
            crux <operation> <options> <arguments>
         where operation is create-index, search, compute-q-values, or
         q-ranker.     
 REVISION: $Revision: 1.2.6.3 $
*/

#include "crux-main.h"

#ifdef CRUX_USE_CUDA
#include "crux_cuda.cu.h"
#endif

#ifdef CRUX_USE_CUDA1
#include "crux_cuda.h"
#endif

#define NUMBER_COMMAND_TYPES 6
static char* command_type_strings[NUMBER_COMMAND_TYPES] =
  {"create-index", "search-for-matches", 
   "compute-q-values", "q-ranker", "percolator", "invalid"};


char* usage_str = "Usage: crux <command> [options] <argument>\n"
"Commands:\n"
"  create-index        Create an index for all peptides in a fasta file.\n"
"  search-for-matches  Search a collection of spectra against a sequence\n"
"                      database, returning a collection of peptide-spectrum\n"
"                      matches (PSMs).\n"
"  compute-q-values    Assign a q-value, which is a statistical confidence\n"
"                      measure that accounts for multiple testing, to each\n"
"                      PSM in a given set.\n" 
"  percolator          Analyze a collection of PSMs to target and decoy\n"
"                      sequences using the percolator algorithm.\n"
/*
"  q-ranker            Analyze a collection of PSMs using the Q-ranker\n"
"                      algorithm.\n"
*/
"Options and arguments:\n"
"  Specific to each command. Type 'crux <command>' to get details.\n"
;
int main(int argc, char** argv){

  // check the syntax for carp <operation>
  if( argc < 2 ){
    fprintf(stderr, usage_str);
    exit(1);
  }

#ifdef CRUX_USE_CUDA
  crux_cuda_initialize();
#endif
#ifdef CRUX_USE_CUDA1
  //initialize cudablas.
  if (!initialize_cudablas())
    fprintf(stderr, "CUDA not initialized, using just CPU\n");
#endif

  // determine the operation
  char* op_string = argv[1];
  COMMAND_T command = string_to_command_type(op_string);

  // call the appropriate function 
  // passing the command line minus the first token ('crux')
  switch(command){
  case INDEX_CMD:
    create_index_main(argc-1, argv+1);
    break;

  case SEARCH_CMD:
    search_main(argc-1, argv+1);
    break;

  case QVALUE_CMD:
    qvalue_main(argc-1, argv+1);
    break;

  case QRANKER_CMD:
    //qranker_main(argc-1, argv+1);
    break;

  case PERCOLATOR_CMD:
    percolator_main(argc-1, argv+1);
    break;

  case INVALID_CMD:
    fprintf(stderr, "Invalid command '%s'\n", op_string);
    fprintf(stderr, usage_str);
    exit(1);
    break;

  }

#ifdef CRUX_USE_CUDA
  crux_cuda_shutdown();
#endif

#ifdef CRUX_USE_CUDA1
  shutdown_cudablas();
#endif

  exit (0);
}// end main

/**
 * \brief Convert a text name of a COMMAND_T to its enumerated type
 * value.  Returns INVALID_CMD for any incorrect string.
 */
COMMAND_T string_to_command_type(char* name){

  int command = convert_enum_type_str( name, 
                                       (int)INVALID_CMD,
                                       command_type_strings, 
                                       NUMBER_COMMAND_TYPES);

  return (COMMAND_T)command;
}















