#include "IonSeries.h"
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

  if (string(peptideB) == string("NULL")) {
    if (posA == -1 || posB == -1) {
      cout<<"Creating linear peptide"<<endl;
      linked_peptide = new LinearPeptide(peptideA);
    } else {
      cout<<"Creating selfloop peptide"<<endl;
      linked_peptide = new SelfLoopPeptide(peptideA, posA-1, posB-1);
    }
  } else {
    cout <<"Creating XLinkPeptide"<<endl;
    linked_peptide = new XLinkPeptide(peptideA, peptideB, posA-1, posB-1);
  }

  cout << "Printing stuff"<<endl;
  cout << "precursor: " << linked_peptide -> getSequenceString() << endl;
  cout << "mass:" << linked_peptide -> getMass(MONO)<<" "<< linked_peptide -> getMass(AVERAGE) <<endl;;
  cout << "charge:" << charge << endl;
  cout << "link mass:"<< linker_mass << endl;
  cout << "print_spectrum:"<< print_spectrum << endl;

  int max_charge = min(get_max_ion_charge_parameter("max-ion-charge"), charge);

  IonConstraint* ion_constraint = new IonConstraint(MONO, max_charge, BY_ION,  false);

  IonSeries* ion_series = new IonSeries(ion_constraint, charge);

  linked_peptide->predictIons(ion_series, charge);

  for (IonIterator ion_iter = ion_series->begin();
    ion_iter != ion_series->end();
    ++ion_iter) {

    Ion* ion = *ion_iter;
    FLOAT_T mz = ion->getMassZ();
    FLOAT_T mass = ion->getMassFromMassZ();
    int charge = ion->getCharge();
    string sequence = linked_peptide->getIonSequence(ion);
    int cleavage_idx = ion->getCleavageIdx();
    cout << mz << "\t"
	 << mass << "\t"
	 << charge << "\t"
         << cleavage_idx << "\t"
	 << sequence << endl;
  }


  delete ion_series;
  delete linked_peptide;

  return 0;
}
