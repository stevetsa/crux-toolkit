#include "ScoreMPSMSpectrum.h"
#include <algorithm>
#include "carp.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "Ion.h"
#include "IonSeries.h"
#include "IonConstraint.h"
#include "SpectrumCollectionFactory.h"

#include "DelimitedFile.h"

#include "MPSM_Scorer.h"

using namespace std;

ScoreMPSMSpectrum::ScoreMPSMSpectrum() {

}

ScoreMPSMSpectrum::~ScoreMPSMSpectrum() {
}


int ScoreMPSMSpectrum::main(int argc, char** argv) {
  /* Define optional and required command line arguments */
  const char* option_list[] = {
    "version",
    "max-ion-charge",
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  const char* argument_list[] = {
    "ms2 file",
    "scan number",
    "peptide sequences",
    "charge states"
  };
  int num_arguments = sizeof(argument_list) / sizeof(char*);

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
  parse_cmd_line_into_params_hash(argc, argv, "score-mpsm-spectrum");

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /* Get Arguments */

  vector<string> peptide_sequences = get_string_vector_parameter("peptide sequences");
  vector<int> charge_states = get_int_vector_parameter("charge states");

  const char* ms2_filename = get_string_parameter("ms2 file");
  int scan_number = get_int_parameter("scan number");

  /* Get Options */
  const char* max_ion_charge = get_string_parameter_pointer("max-ion-charge");

  // check peptide sequence
  for (int idx=0;idx<peptide_sequences.size();idx++) {
    if(!valid_peptide_sequence(peptide_sequences[idx].c_str())){
      carp(CARP_FATAL, "The peptide sequence '%s' is not valid", 
           peptide_sequences[idx].c_str());
    }
  }

  int max_charge = charge_states[0];
  for (unsigned int idx=1;idx < charge_states.size();idx++) {
    max_charge = max(charge_states[idx], max_charge);
  }

  max_charge = min(max_charge, get_max_ion_charge_parameter("max-ion-charge"));

  //get spectrum
  
  SpectrumCollection* collection = 
    SpectrumCollectionFactory::create(ms2_filename);

  Spectrum* spectrum = collection->getSpectrum(scan_number);

  if (spectrum == NULL) {
    carp(CARP_FATAL, "Could not find scan number %i", scan_number);
  }


  // create ion_constraint
  MASS_TYPE_T frag_masses = get_mass_type_parameter("fragment-mass");
  
  IonConstraint* ion_constraint = 
    IonConstraint::newIonConstraintSmart(XCORR, max_charge);

  IonSeries* ion_series = new IonSeries(ion_constraint, max_charge);

  MPSM_Scorer::createIonSeries(peptide_sequences, charge_states, ion_series);

  Scorer* scorer = new Scorer(XCORR);
  
  FLOAT_T xcorr = scorer->scoreSpectrumVIonSeries(spectrum, ion_series);

  delete scorer;
  delete ion_series;
  IonConstraint::free(ion_constraint);

  ion_constraint = 
    IonConstraint::newIonConstraintSmart(SP, max_charge);

  ion_series = new IonSeries(ion_constraint, max_charge);

  MPSM_Scorer::createIonSeries(peptide_sequences, charge_states, ion_series);

  scorer = new Scorer(SP);
  FLOAT_T sp = scorer->scoreSpectrumVIonSeries(spectrum, ion_series);

  delete scorer;
  delete ion_series;
  IonConstraint::free(ion_constraint);

  // print settings
/*  
  for (int idx=0;idx < peptide_sequences.size();idx++) {
    printf("# INDEX: %d\n",idx);  
    printf("# PEPTIDE: %s\n",peptide_sequences[idx].c_str());
    printf("# AVERAGE: %f MONO:%f\n",
      calc_sequence_mass(peptide_sequences[idx].c_str(), AVERAGE),
      calc_sequence_mass(peptide_sequences[idx].c_str(), MONO));
    printf("# CHARGE: %d\n", charge_states[idx]);
  }
  printf("# MAX-ION-CHARGE: %s\n", max_ion_charge);
*/
  printf("%f\n", xcorr);
  printf("%f\n", sp);

  carp(CARP_INFO, "score-mpsm-spectrum finished");
  return 0;
}

string ScoreMPSMSpectrum::getName() {
  return "score-mpsm-spectrum";
}

string ScoreMPSMSpectrum::getDescription() {
  return 
    "Given a list a peptides, charges, calculates the xcorr and sp for a scan in an ms2file";
}
