/**
 * \file CruxHK2MS2Application.cpp 
 * \brief Given a ms1 and ms2 file, run hardklor followed by the bullseye algorithm.
 *****************************************************************************/
#include "CruxHK2MS2Application.h"
#include "CruxHardklorApplication.h"
#include "DelimitedFileWriter.h"

#include "crux-utils.h"

using namespace std;

/**
 * \returns a blank CruxHK2MS2Application object
 */
CruxHK2MS2Application::CruxHK2MS2Application() {

}

/**
 * Destructor
 */
CruxHK2MS2Application::~CruxHK2MS2Application() {
}

/**
 * main method for CruxHK2MS2Application
 */
int CruxHK2MS2Application::main(int argc, char** argv) {

   /* Define optional command line arguments */
  const char* option_list[] = {
    "fileroot",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity",
    "peak-type"
  };

  int num_options = sizeof(option_list) / sizeof(char*);

  /* Define required command line arguments */
  const char* argument_list[] = {"MS2 spectra"};
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  /* Initialize the application */

  initialize(argument_list, num_arguments,
    option_list, num_options, argc, argv);

  /* Get parameters. */
  string hardklor_output;
  string input_ms2 = get_string_parameter("MS2 spectra");
  
  string output_stem = make_file_path("hk2ms2");

  bool overwrite = get_boolean_parameter("overwrite");

  hardklor_output = string(get_string_parameter_pointer("hardklor-file"));

  if (hardklor_output == "__NULL_STR") {
    hardklor_output = make_file_path("hardklor.mono.txt");
    if ((overwrite) || (!file_exists(hardklor_output))) {
      carp(CARP_DEBUG,"Calling hardklor");
      bool ret = CruxHardklorApplication::main(input_ms2);
      if (ret != 0) {
        carp(CARP_ERROR, "Hardklor failed:%d", ret);
        return ret;
      }
    }
  }

  DEISOTOPE_PEAKS_T deisotope_output_peak = get_deisotope_peaks_parameter("peak-type");


  /* build argument list */
  vector<string> hk2ms2_args_vec;
  hk2ms2_args_vec.push_back("hk2ms2");
    
  hk2ms2_args_vec.push_back(hardklor_output);
  hk2ms2_args_vec.push_back(input_ms2);
  hk2ms2_args_vec.push_back(output_stem);

  /* add flags */
  switch (deisotope_output_peak) {
    case MONOISOTOPIC_DEISOTOPE_PEAKS:
      //default for hk2ms2, no flag needed.
      break;
    case MASS_PLUS_H_DEISOTOPE_PEAKS:
      hk2ms2_args_vec.push_back("-h");
      break;
    case MASS_TO_CHARGE_DEISOTOPE_PEAKS:
      hk2ms2_args_vec.push_back("-z");
      break;

    case INVALID_DEISOTOPE_PEAKS:
    case NUMBER_DEISOTOPE_PEAKS:
      carp(CARP_FATAL, "Invalid deisotope peaks type!");
      
  }

  /* build argv line */
  int hk2ms2_argc = hk2ms2_args_vec.size();

  char** hk2ms2_argv = new char*[hk2ms2_argc];

  hk2ms2_argv[0] = (char*)hk2ms2_args_vec[0].c_str();
  for (int idx = 1;idx < hk2ms2_argc ; idx++) {
    hk2ms2_argv[idx] = (char*)hk2ms2_args_vec[idx].c_str();
    carp(CARP_DEBUG, "hk2ms2_argv[%d]=%s", idx, hk2ms2_argv[idx]);
  }

  /* Call hk2ms2Main */
  int ret = hk2ms2Main(hk2ms2_argc, hk2ms2_argv);

  delete []hk2ms2_argv;

  return ret;
}

/**
 * \returns the command name for CruxHK2MS2Application
 */
string CruxHK2MS2Application::getName() {
  return "deisotope-ms2";
}

/**
 * \returns the description for CruxHK2MS2Application
 */
string CruxHK2MS2Application::getDescription() {

  return "Runs centroiding/deisotoping algorithm using hardklor";
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool CruxHK2MS2Application::needsOutputDirectory() {
  return true;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
