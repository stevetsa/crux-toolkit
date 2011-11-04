/**
 * \file CruxHardklorApplication.cpp 
 * \brief Given a delimited file and a column-name, print out statistics
 * for that column (n, min, max, sum, average, stddev, median).
 *****************************************************************************/
#include "CruxHardklorApplication.h"

#include "DelimitedFile.h"

using namespace std;


/**
 * \returns a blank CruxHardklorApplication object
 */
CruxHardklorApplication::CruxHardklorApplication() {

}

/**
 * Destructor
 */
CruxHardklorApplication::~CruxHardklorApplication() {
}

/**
 * main method for CruxHardklorApplication
 */
int CruxHardklorApplication::main(int argc, char** argv) {

   /* Define optional command line arguments */
  const char* option_list[] = {
    "delimiter",
    "header",
    "precision",
    "verbosity"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"tsv file", "column name"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */
  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

  carp(CARP_FATAL, "Not implemented yet!");



  return 0;
}

/**
 * \returns the command name for CruxHardklorApplication
 */
string CruxHardklorApplication::getName() {
  return "hardklor";
}

/**
 * \returns the description for CruxHardklorApplication
 */
string CruxHardklorApplication::getDescription() {

  return "Runs hardklor";
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
