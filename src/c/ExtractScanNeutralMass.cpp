#include "ExtractScanNeutralMass.h"
#include "SpectrumCollectionFactory.h"

#include <iostream>

//#include "print-processed-spectra.h"

using namespace std;

ExtractScanNeutralMass::ExtractScanNeutralMass() {

}

ExtractScanNeutralMass::~ExtractScanNeutralMass() {
}


int ExtractScanNeutralMass::main(int argc, char** argv) {

  const char* argument_list[] = { "ms2 file" };

  int num_arguments = sizeof(argument_list) / sizeof(char*);

  const char* option_list[] = {
    "verbosity"
  };

  int num_options = sizeof(option_list) / sizeof(char*);

      // Verbosity level for set-up/command line reading 
  set_verbosity_level(CARP_WARNING);

  // Initialize parameter.c and set default values
  initialize_parameters();

  // Define optional and required arguments 
  select_cmd_line_options(option_list, num_options);
  select_cmd_line_arguments(argument_list, num_arguments);
  
  // Parse the command line, including optional params file
  // Includes syntax, type, and bounds checking, dies on error 
  const char* cmd_name = this->getName().c_str();
  char* full_cmd = cat_string("crux ", cmd_name);
  parse_cmd_line_into_params_hash(argc, argv, cmd_name);
  free(full_cmd);


  const char* ms2_file = get_string_parameter_pointer("ms2 file");

  SpectrumCollection* spectra = SpectrumCollectionFactory::create(ms2_file);

  spectra->parse();

  cout << "scan\tcharge\tneutral mass" << endl; 

  FilteredSpectrumChargeIterator* spectrum_iterator =
    new FilteredSpectrumChargeIterator(spectra);

  if (spectrum_iterator == NULL) {
    carp(CARP_FATAL, "Could not create spectrum iterator");
  }

  while (spectrum_iterator->hasNext()) {

    SpectrumZState cur_zstate;
    
    Spectrum* cur_spectrum = 
      spectrum_iterator->next(cur_zstate);

    cout << cur_spectrum->getFirstScan() << "\t"
         << cur_zstate.getCharge() << "\t"
         << cur_zstate.getNeutralMass() << "\t" 
         << cur_zstate.getRTime() << endl;

  }

  
  delete spectra;
  
  


  return 0;
}

string ExtractScanNeutralMass::getName() {
  return "extract-scan-neutral-mass";
}

string ExtractScanNeutralMass::getDescription() {
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

