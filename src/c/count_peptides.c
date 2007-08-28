/*****************************************************************************
 * \file generate_peptides
 * AUTHOR: Chris Park
 * CREATE DATE: July 17 2006
 * DESCRIPTION: Given a protein fasta sequence database as input, generate a list of peptides in 
 *              the database that meet certain criteria (e.g. mass, length, trypticity) as output.
 * REVISION: 
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include "carp.h"
#include "utils.h"
#include "crux-utils.h"
#include "objects.h"
#include "mass.h"
#include "peptide.h"
#include "peptide_src.h"
#include "protein.h"
#include "database.h"
#include "parse_arguments.h"
#include "parameter.h"
#include "index.h"
#include "generate_peptides_iterator.h"

/**
 * when wrong command is seen carp, and exit
 */
void wrong_command(char* arg){
  char* usage = parse_arguments_get_usage("generate_peptides");
  carp(CARP_FATAL, "incorrect argument %s", arg);
  fprintf(stderr, "%s", usage);
  free(usage);
  exit(1);
}

/**
 * progress indicator, displays a spinning | to standout
 */
void show_progress(int* num){
  putc('\b', stderr);
  if(*num / 150 == 37){
    putc('|', stderr);
    fflush(stderr);
  }
  else if(*num / 150 == 74){
    putc('/', stderr);
    fflush(stderr);
  }
  else if(*num / 150 == 111){
    putc('-', stderr);
    fflush(stderr);
  }
  else if(*num / 150 == 149){
    putc('\\', stderr);
    fflush(stderr);
    *num = 0;
    return;
  }
  ++*num;
}

int main(int argc, char** argv){

  /* Set default values for any options here */
  int flag_opt = FALSE;
  int trypticity_opt = FALSE;
  double min_mass = 200;
  double max_mass = 7200;
  int min_length = 6;
  int max_length = 50;
  char* cleavages = "tryptic"; 
  char* isotopic_mass = "average" ;
  int  verbosity = CARP_FATAL;
  char* redundancy = "redundant";
  char* use_index = "F";
  char* parameter_file = "crux_parameter";

  int missed_cleavages = FALSE;
  char* sort = "none";      // mass, length, lexical, none  
  char * in_file = NULL;
  char * error_message;
  int result = 0;

  /* Define optional command line arguments */ 
  
  parse_arguments_set_opt(
    "min-mass", 
    "The minimum neutral mass of the peptides to output.", 
    (void *) &min_mass, 
    DOUBLE_ARG);

  parse_arguments_set_opt(
    "max-mass", 
    "The maximum neutral mass of the peptides to output.", 
    (void *) &max_mass, 
    DOUBLE_ARG);

  parse_arguments_set_opt(
    "min-length", 
    "The minimum length of the peptides to output.",
    (void *) &min_length, 
    INT_ARG);

  parse_arguments_set_opt(
    "max-length", 
    "The maximum length of the peptides to output. maximum limit = 255.",
    (void *) &max_length, 
    INT_ARG);

  parse_arguments_set_opt(
    "cleavages", 
    "Type of cleavages to allow. tryptic|partial|all.", 
    (void *) &cleavages, 
    STRING_ARG);

  parse_arguments_set_opt(
    "missed-cleavages", 
    "Allow missed cleavage sites with in a peptide. ",
    (void *) &missed_cleavages, 
    FLAG_ARG);
  
  parse_arguments_set_opt(
    "sort", 
    "Specify the order in which peptides are printed to standard output. none|mass|length|lexical.", 
    (void *) &sort, 
    STRING_ARG);

  parse_arguments_set_opt(
    "isotopic-mass", 
    "Specify the type of isotopic masses to use when calculating the peptide mass. average|mono.",
    (void *) &isotopic_mass, 
    STRING_ARG);

  parse_arguments_set_opt(
    "verbosity", 
    "Specify the verbosity of the current processes from 0-100.",
    (void *) &verbosity, 
    INT_ARG);

  parse_arguments_set_opt(
    "redundancy", 
    "Specify whether peptides that come from different proteins yet with identical sequences should appear on separate lines or on the same line. redundant|unique.",
    (void *) &redundancy, 
    STRING_ARG);
  
  parse_arguments_set_opt(
    "use-index", 
    "Specify whether a pre-computed on-disk index should be used for retrieving the peptides. T|F",
    (void *) &use_index, 
    STRING_ARG);

  parse_arguments_set_opt(
    "parameter-file",
    "The crux parameter file to parse parameter from.",
    (void *) &parameter_file,
    STRING_ARG);

  parse_arguments_set_opt(
    "output-sequence",
    "Output the peptide sequence as well as the protein id and start and stop.",
    (void *) &flag_opt,
    FLAG_ARG);

  parse_arguments_set_opt(
    "output-trypticity",
    "Output the peptide trypticiy in output.",
    (void *) &trypticity_opt,
    FLAG_ARG);

  /* Define required command line arguments */
  parse_arguments_set_req(
    "fasta-file",
    "The name of the file (in fasta format) from which to parse proteins.",
    (void *) &in_file, STRING_ARG);


  /* Parse the command line */
  if (parse_arguments(argc, argv, 0)) {
    GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator = NULL;

    //set verbosity
    if(CARP_FATAL <= verbosity && verbosity <= CARP_MAX){
      set_verbosity_level(verbosity);
    }
    else{
      wrong_command("verbosity");
    }

    //first, parse paramter file
    //Next, updates the parameter files with command line options
    parse_update_parameters(parameter_file);

    //parameters are now confirmed, can't be changed
    parameters_confirmed();

    //create peptide interator
    peptide_iterator = new_generate_peptides_iterator();

    print_peptide_count(peptide_iterator);

    free_generate_peptides_iterator(peptide_iterator);

    free_parameters();
    exit(0);
  }
  else {
    char* usage = parse_arguments_get_usage("generate_peptides");
    result = parse_arguments_get_error(&error_message);
    fprintf(stderr, "Error in command line. Error # %d\n", result);
    fprintf(stderr, "%s\n", error_message);
    fprintf(stderr, "%s", usage);
    free(usage);
  }
  exit(0);
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
