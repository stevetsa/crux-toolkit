#include "PrintProcessedSpectra.h"

using namespace std;

PrintProcessedSpectra::PrintProcessedSpectra() {

}

PrintProcessedSpectra::~PrintProcessedSpectra() {
}


int PrintProcessedSpectra::main(int argc, char** argv) {

  // Define optional command line arguments
  const char* option_list[] = { 
    "verbosity",
    "parameter-file", 
    "overwrite"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  // Define required command line arguments
  const char* argument_list[] = { "ms2 file", 
                                  "output file"};

  int num_arguments = sizeof(argument_list) / sizeof(char*);
  // For output of parameter parsing
  set_verbosity_level(CARP_ERROR);  

  // set up parameters and their defaults in parameter.c
  initialize_parameters();

  // Define optional and required command line arguments
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments);

  // Parse the command line, including the optional params file
  parse_cmd_line_into_params_hash(argc, argv, "crux print-processed-spectra");

  // Get arguments and options
  const char* input_ms2_name  = get_string_parameter_pointer("ms2 file");
  char* output_ms2_name = get_string_parameter("output file");
  prefix_fileroot_to_name(&output_ms2_name);
  const char* output_dir = get_string_parameter_pointer("output-dir");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");

  // open output file
  create_output_directory(output_dir, overwrite);
  FILE* output_ms2 = create_file_in_path(output_ms2_name,
                                         output_dir,
                                         overwrite);
  // open input file
  SpectrumCollection* spectra = new SpectrumCollection(input_ms2_name);
  if( spectra == NULL ){
    carp(CARP_FATAL, "Could not read spectra from %s.", input_ms2_name);
  }

  spectra->parse();
  carp(CARP_DEBUG, "Found %d spectra in file.", 
       spectra->getNumSpectra());

  // write header to output file
  char* header = spectra->getComment();
  fprintf(output_ms2, "%s", header);
  fprintf(output_ms2, "H\tComment\tSpectra processed as for Xcorr\n");

  // create iterator for getting spectra
  FilteredSpectrumChargeIterator* spectrum_iterator =
    new FilteredSpectrumChargeIterator(spectra);

  if( spectrum_iterator == NULL ){
    carp(CARP_FATAL, "Could create spectrum iterator");
  }

  // loop over all spectra, process, print
  while(spectrum_iterator->hasNext()){
    int cur_charge = 0;
    Spectrum* cur_spectrum = 
      spectrum_iterator->next(&cur_charge);

    carp(CARP_DETAILED_INFO, "Processing spectrum %d charge %d.",
         cur_spectrum->getFirstScan(), cur_charge);

    // change the peak values
    FLOAT_T* intensities = NULL;
    int max_mz_bin = 0;
    get_processed_peaks(cur_spectrum, cur_charge,
                        &intensities, &max_mz_bin);

    // print processed spectrum
    cur_spectrum->printProcessedPeaks(cur_charge, 
                                        intensities, max_mz_bin,
                                        output_ms2);
  }

  // close output file
  delete spectra;
  fclose(output_ms2);

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux print-processed-spectra.");

  return(0);
}

string PrintProcessedSpectra::getName() {
  return "print-processed-spectra";
}

string PrintProcessedSpectra::getDescription() {
  return 
    "Write a new ms2 file with all of the same spectra "
    "with only the peaks used for computing xcorr.";
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

