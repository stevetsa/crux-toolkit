#include "PredictMPSMIons.h"
#include <algorithm>
#include "carp.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "Ion.h"
#include "IonSeries.h"
#include "IonConstraint.h"

#include "DelimitedFile.h"

#include "MPSM_Scorer.h"

using namespace std;

PredictMPSMIons::PredictMPSMIons() {

}

PredictMPSMIons::~PredictMPSMIons() {
}


int PredictMPSMIons::main(int argc, char** argv) {
  /* Define optional and required command line arguments */
  const char* option_list[] = {
    "version",
    "primary-ions",
    "precursor-ions",
    "neutral-losses",
    "isotope",
    "flanking",
    "max-ion-charge",
    "nh3",
    "h2o"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  const char* argument_list[] = {
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
  parse_cmd_line_into_params_hash(argc, argv, "predict-mpsm-ions");

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /* Get Arguments */

  vector<string> peptide_sequences = get_string_vector_parameter("peptide sequences");
  vector<int> charge_states = get_int_vector_parameter("charge states");

  /* Get Options */
  ION_TYPE_T ion_type = get_ion_type_parameter("primary-ions");
  BOOLEAN_T use_precursor_ions = get_boolean_parameter("precursor-ions");
  int isotope_count = get_int_parameter("isotope");
  BOOLEAN_T is_flanking = get_boolean_parameter("flanking");
  const char* max_ion_charge = get_string_parameter_pointer("max-ion-charge");
  int nh3_count = get_int_parameter("nh3");
  int h2o_count = get_int_parameter("h2o");

  int neutral_loss_count[MAX_MODIFICATIONS];
  BOOLEAN_T is_modification = FALSE;

  // check peptide sequence
  for (int idx=0;idx<peptide_sequences.size();idx++) {
    if(!valid_peptide_sequence(peptide_sequences[idx].c_str())){
      carp(CARP_FATAL, "The peptide sequence '%s' is not valid", 
           peptide_sequences[idx].c_str());
    }
  }

   // neutral_losses
   // initialize
   int modification_idx = 0;
   for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
     neutral_loss_count[modification_idx] = 0;
   }

   is_modification = (nh3_count || 
                      h2o_count || 
                      is_flanking);

   neutral_loss_count[NH3] = nh3_count;
   neutral_loss_count[H2O] = h2o_count;
   neutral_loss_count[FLANK] = (int)is_flanking;
   neutral_loss_count[ISOTOPE] = isotope_count;

  

  int max_charge = charge_states[0];
  for (unsigned int idx=1;idx < charge_states.size();idx++) {
    max_charge = max(charge_states[idx], max_charge);
  }

  max_charge = min(max_charge, get_max_ion_charge_parameter("max-ion-charge"));
  // create ion_constraint
  MASS_TYPE_T frag_masses = get_mass_type_parameter("fragment-mass");
  
  IonConstraint* ion_constraint = 
  //  new_ion_constraint(MONO, max_charge, ion_type, use_precursor_ions);
    new IonConstraint(frag_masses, max_charge, ion_type, use_precursor_ions);

   
   // set ion_constraint3 modification counts, if modifications should occur
  if(is_modification){
    ion_constraint->setModification(NH3, neutral_loss_count[NH3]);
    ion_constraint->setModification(H2O, neutral_loss_count[H2O]);
    ion_constraint->setModification(ISOTOPE, neutral_loss_count[ISOTOPE]);
    ion_constraint->setModification(FLANK, neutral_loss_count[FLANK]);
  }

  IonSeries* ion_series = new IonSeries(ion_constraint, max_charge);

  MPSM_Scorer::createIonSeries(peptide_sequences, charge_states, ion_series);

  // print settings
  
  for (int idx=0;idx < peptide_sequences.size();idx++) {
    printf("# INDEX: %d\n",idx);  
    printf("# PEPTIDE: %s\n",peptide_sequences[idx].c_str());
    printf("# AVERAGE: %f MONO:%f\n",
      calc_sequence_mass(peptide_sequences[idx].c_str(), AVERAGE),
      calc_sequence_mass(peptide_sequences[idx].c_str(), MONO));
    printf("# CHARGE: %d\n", charge_states[idx]);
  }
  printf("# MAX-ION-CHRAGE: %s\n", max_ion_charge);
  printf("# NH3 modification: %d\n", neutral_loss_count[NH3]);
  printf("# H2O modification: %d\n", neutral_loss_count[H2O] );
  printf("# ISOTOPE modification: %d\n", neutral_loss_count[ISOTOPE] );
  printf("# FLANK modification: %d\n", neutral_loss_count[FLANK]);
  // print ions
  ion_series->print(stdout);

  // free
  IonConstraint::free(ion_constraint);

  delete ion_series;
  carp(CARP_INFO, "predict-mpsm-ions finished");
  return 0;
}

string PredictMPSMIons::getName() {
  return "predict-mpsm-ions";
}

string PredictMPSMIons::getDescription() {
  return 
    "Given a list a peptides and charges, predicts the ions ";
}
