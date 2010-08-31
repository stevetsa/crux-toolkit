#include "ion_series.h"
#include "MatchCandidate.h"
#include "XLinkPeptide.h"
#include "LinearPeptide.h"
#include "SelfLoopPeptide.h"

#include <iostream>
#include <fstream>

#define NUM_ARGUMENTS 6
#define NUM_OPTIONS 4
using namespace std;

int main(int argc, char** argv) {

  char* peptideA = NULL;
  char* peptideB = NULL;
  int posA = 0;
  int posB = 0;
  FLOAT_T linker_mass = 0;
  int charge = 1; 
  BOOLEAN_T print_spectrum = FALSE;


  /* Verbosity level for set-up/command line reading */
  set_verbosity_level(CARP_ERROR);
  
  /* Define optional command line arguments */
  int num_options = NUM_OPTIONS;
  const char* option_list[NUM_OPTIONS] = {
    "verbosity",
    "version",
    "xcorr-use-flanks",
    "print-theoretical-spectrum"
  };

  /* Define required command line arguments */
  int num_arguments = NUM_ARGUMENTS;
  const char* argument_list[NUM_ARGUMENTS] = {"peptide A",
                                              "peptide B",
					      "pos A",
					      "pos B",
					      "charge state",
					      "link mass"};



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
  parse_cmd_line_into_params_hash(argc, argv, "xlink-predict-peptide-ions");

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /* Get Arguments */
  linker_mass = get_double_parameter("link mass");
  charge = get_int_parameter("charge state");

  peptideA = get_string_parameter("peptide A");
  peptideB = get_string_parameter("peptide B");
  
  posA = get_int_parameter("pos A");
  posB = get_int_parameter("pos B");

  print_spectrum = get_boolean_parameter("print-theoretical-spectrum");

  XLinkPeptide::setLinkerMass(linker_mass);

  MatchCandidate* linked_peptide = NULL;

  cout <<"creating peptide"<<endl;

  if (string(peptideB) == "NULL") {
    if (posA == -1 || posB == -1) {
      linked_peptide = new LinearPeptide(peptideA);
    } else {
      linked_peptide = new SelfLoopPeptide(peptideA, posA-1, posB-1);
    }
  } else {
    cout <<"Creating XLinkPeptide"<<endl;
    linked_peptide = new XLinkPeptide(peptideA, peptideB, posA-1, posB-1);
  }

  cout << "Printing stuff"<<endl;
  cout << "precursor: " << linked_peptide -> getSequenceString() << endl;
  cout << "mass:" << linked_peptide -> getMass()<<endl;;
  cout << "charge:" << charge << endl;
  cout << "link mass:"<< linker_mass << endl;
  cout << "print_spectrum:"<< print_spectrum << endl;

  int max_charge = min(get_max_ion_charge_parameter("max-ion-charge"), charge);

  ION_CONSTRAINT_T* ion_constraint = new_ion_constraint(MONO, max_charge, BY_ION,  FALSE);

  ION_SERIES_T* ion_series = new_ion_series_generic(ion_constraint, charge);

  linked_peptide->predictIons(ion_series, charge);

  ION_ITERATOR_T* ion_iter = new_ion_iterator(ion_series);

  while (ion_iterator_has_next(ion_iter)) {
    ION_T* ion = ion_iterator_next(ion_iter);
    FLOAT_T mz = get_ion_mass_z(ion);
    FLOAT_T mass = get_ion_mass_from_mass_z(ion);
    int charge = get_ion_charge(ion);
    string sequence = linked_peptide->getIonSequence(ion);
    int cleavage_idx = get_ion_cleavage_idx(ion);
    cout << mz << "\t"
	 << mass << "\t"
	 << charge << "\t"
         << cleavage_idx << "\t"
	 << sequence << endl;
    
    

  }

  free_ion_iterator(ion_iter);
  free_ion_series(ion_series);

  delete linked_peptide;

  return 0;
}
