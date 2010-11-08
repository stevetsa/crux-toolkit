/*************************************************************************//**
 * \file predict_peptide_ions.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 10/05 2006
 * \brief Object for given a peptide and a charge state, predict
 * the ions 
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include "carp.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"

#include "AKlammerRetentionPredictor.h"
#include "DelimitedFileReader.h"



int main(int argc, char** argv){

  /* Define optional and required command line arguments */

  const char* option_list[] = {
    "version",
    "aklammer-rtime-nterm",
    "aklammer-rtime-cterm",
    "aklammer-rtime-tryptic"
  };
  int num_options = sizeof(option_list)/sizeof(const char*);
 
  const char* argument_list[] = {
    "result file",
  };

  int num_arguments = sizeof(argument_list)/sizeof(const char*);
  

  /* for debugging of parameter processing */
  //set_verbosity_level( CARP_DETAILED_DEBUG );
  set_verbosity_level( CARP_ERROR );

  /* Set default values for parameters in parameter.c */
  initialize_parameters();

  /* Define optional and required command line arguments */
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments);

  /* Parse the command line, including the optional params file */
  /* does sytnax, type, bounds checking and dies if neccessessary */
  parse_cmd_line_into_params_hash(argc, argv, "crux-predict-peptide-ions");

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /* Get Arguments */

  const char* result_file_name = get_string_parameter_pointer("result file");

  /* Get Options */

  /* Train */

  DelimitedFileReader result_file(result_file_name, true);

  AKlammerRetentionPredictor::train(result_file);

   carp(CARP_INFO, "crux-train-aklammer finished");
 exit(0);
}
