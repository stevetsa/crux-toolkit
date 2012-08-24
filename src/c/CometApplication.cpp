/**
 * \file CometApplication.cpp 
 * \brief Runs hardklor
 *****************************************************************************/
#include "CometApplication.h"
#include "DelimitedFileWriter.h"
#include "DelimitedFile.h"
#include "Common.h"

using namespace std;

/**
 * \returns a blank CometApplication object
 */
CometApplication::CometApplication() {

}

/**
 * Destructor
 */
CometApplication::~CometApplication() {
}

/**
 * main method for CometApplication
 */
int CometApplication::main(int argc, char** argv) {

   /* Define optional command line arguments */
  const char* option_list[] = {
    "fileroot",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity"
  };

  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */

  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

  

  return comet_main(argc, argv);

}

/**
 * \returns the command name for CometApplication
 */
string CometApplication::getName() {
  return "comet";
}

/**
 * \returns the description for CometApplication
 */
string CometApplication::getDescription() {

  return "Runs comet";
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool CometApplication::needsOutputDirectory() {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
