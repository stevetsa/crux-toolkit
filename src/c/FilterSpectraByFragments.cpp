/**
 * \file PrintProcessedSpectra.cpp
 *
 * AUTHOR: Barbara Frewen
 * CREATE DATE: September 18, 2009
 *
 * DESCRIPTION: Main method for the print-processed-spectra command.
 *              For every spectrum in an ms2 file, process as for
 *              xcorr and print peaks in ms2 format to new file.
 */

#include "FilterSpectraByFragments.h"
#include "SpectrumCollection.h"
#include "SpectrumCollectionFactory.h"

using namespace std;

/**
 * \returns a blank FilterSpectraByFragments object
 */
FilterSpectraByFragments::FilterSpectraByFragments() {

}

/**
 * Destructor
 */
FilterSpectraByFragments::~FilterSpectraByFragments() {
}

bool FilterSpectraByFragments::hasFragments(
  Spectrum* spectrum, 
  vector<double>& mz_fragment_list,
  double mz_tolerance) {

  //try to find each peak in the spectrum.
  for (unsigned int frag_idx = 0;
    frag_idx < mz_fragment_list.size();
    ++frag_idx) {
    
    PEAK_T* peak = spectrum->getNearestPeak(mz_fragment_list[frag_idx], mz_tolerance);
    if (peak == NULL) {
      return false;
    }
  }

  return true;

}



/**
 * main method for FilterSpectraByFragments
 */
int FilterSpectraByFragments::main(int argc, char** argv) {
  
  // Define optional command line arguments
  const char* option_list[] = { 
    "verbosity",
    "parameter-file", 
    "overwrite"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  // Define required command line arguments
  const char* argument_list[] = { "ms2 file", "fragment-masses", "mz-bin-width"}; 
  int num_arguments = sizeof(argument_list) / sizeof(char*);
  
  initialize(argument_list, 
             num_arguments,
             option_list, 
             num_options,
             argc, argv);

  // Get arguments and options
  const char* input_name  = get_string_parameter_pointer("ms2 file");
  char* output_ms2_name = get_string_parameter("output file");
  prefix_fileroot_to_name(&output_ms2_name);

  double mz_tolerance = get_double_parameter("mz-bin-width");
  vector<double> mz_fragment_list = get_double_vector_parameter("fragment-masses");


  // open input file
  SpectrumCollection* spectra = SpectrumCollectionFactory::create(input_name);
  if( spectra == NULL ){
    carp(CARP_FATAL, "Could not read spectra from %s.", input_name);
  }

  spectra->parse();
  carp(CARP_DEBUG, "Found %d spectra in file.", 
       spectra->getNumSpectra());

  // write header to output file
  //const char* header = spectra->getComment();

  
  //fprintf(stdout, "%s", header);
  fprintf(stdout, "H\tComment\tSpectra filtered by fragment m/z\n");

  // loop over all spectra, process, print

  int kept_spectra = 0;
  int total_spectra = spectra->getNumSpectra();

  for (SpectrumIterator spectrum_iter = spectra->begin();
    spectrum_iter != spectra->end();
    ++spectrum_iter) {

    Spectrum* cur_spectrum = *spectrum_iter;

    if (hasFragments(cur_spectrum, mz_fragment_list, mz_tolerance)) {
      cur_spectrum->print(stdout);
      kept_spectra++;
    }
  }
  delete spectra;

  carp(CARP_INFO,"Kept %i out of %i spectra", kept_spectra, total_spectra);

  return(0);

}

/**
 * \returns the command name for FilterSpectraByFragments
 */
string FilterSpectraByFragments::getName() {
  return "filter-spectra-by-fragments";
}

/**
 * \returns the description for FilterSpectraByFragments
 */
string FilterSpectraByFragments::getDescription() {
  return 
    "Print out a ms2 file that have spectra that contain "
    "certain m/z fragments ";
}

/**
 * \returns the file stem of the application, default getName.
 */
string FilterSpectraByFragments::getFileStem() {
  return "filter-spectra-by-fragments";
}

COMMAND_T FilterSpectraByFragments::getCommand() {
  return MISC_COMMAND;
}

bool FilterSpectraByFragments::needsOutputDirectory() {
  return false;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

