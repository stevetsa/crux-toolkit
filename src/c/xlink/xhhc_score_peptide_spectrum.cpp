#include "xhhc.h"
#include "xhhc_ion_series.h"
#include "xhhc_scorer.h"
#include "objects.h"
#include <math.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#define bin_width_mono 1.0005079


//int noflanks = 0;
//int normalize = 0;

int main(int argc, char** argv){

  // required variables
  char* ms2_file = NULL;
  int scan_num = 0;
 
  //char* positions = NULL;
  char* sequenceA = NULL;
  char* sequenceB = NULL;
  int positionA = 0;
  int positionB = 0;

  // optional variables
  char* linker_mass_string = "0";
  int print_spectrums_int = 0;
  int charge = 2;
  //char* type = "xcorr";
  char* parameter_file = "crux.params";
  int  verbosity = CARP_ERROR;
  // parsing variables
  int result = 0;
  char * error_message; 
  /* Define optional command line arguments */ 

  parse_arguments_set_opt(
    "charge", 
    "The peptide charge. 1|2|3",
    (void *) &charge, 
    INT_ARG);

  parse_arguments_set_opt(
    "parameter-file",
    "The crux parameter file to parse parameter from.",
    (void *) &parameter_file,
    STRING_ARG);


  parse_arguments_set_opt(
    "linker-mass", 
    "linker mass and linker modifications. Default 0.",
    (void *) &linker_mass_string, 
    STRING_ARG);
/*
  parse_arguments_set_opt(
    "no-flanks", 
    "graph theoretical spectrum without flanks",
    (void *) &noflanks, 
    INT_ARG);
  
  parse_arguments_set_opt(
    "normalize", 
    "normalize theoretical spectrum intensities",
    (void *) &normalize, 
    INT_ARG);
*/
  parse_arguments_set_opt(
    "print-spectrums", 
    "",
    (void *) &print_spectrums_int, 
    INT_ARG);
  parse_arguments_set_opt(
    "verbosity", 
    "Specify the verbosity of the current processes from 0-100.",
    (void *) &verbosity, 
    INT_ARG);
  
 /* Define required command line arguments */
  parse_arguments_set_req(
    "peptide-sequence-alpha", 
    "The first peptide sequence.", 
    (void *) &sequenceA, 
    STRING_ARG);

  parse_arguments_set_req(
    "peptide-sequence-beta", 
    "The second peptide sequence.", 
    (void *) &sequenceB, 
    STRING_ARG);

 parse_arguments_set_req(
	"position A",
 	"zero based index of linker on peptide A [0, length(A)-1]",
	(void *) &positionA,
	INT_ARG);

 parse_arguments_set_req(
	"position B",
 	"zero based index of linker on peptide B [0, length(B)-1]",
	(void *) &positionB,
	INT_ARG);
/*
  parse_arguments_set_req(
    "link-at", 
    "Specific link positions on peptide A and B ex. 4,0",
    (void *) &positions, 
    STRING_ARG);
*/  
  parse_arguments_set_req(
    "scan-number", 
    "The scan number for the MS-MS spectrum to extract from the ms2 file. This is an integer in the range [1, 100000], and uniquely identifies a particular MS-MS spectrum within an .ms2 file.",
    (void *) &scan_num, INT_ARG);

  parse_arguments_set_req(
    "ms2-filename", 
    "A file containing multiple MS-MS spectra in .ms2 format.",
    (void *) &ms2_file,
    STRING_ARG);

  /* Parse the command line */
  if (parse_arguments(argc, argv, 0)) {
    // parsed arguments
    bool print_spectrums = (bool) print_spectrums_int;    
    SPECTRUM_T* spectrum = NULL;
    SPECTRUM_COLLECTION_T* collection = NULL;
    FLOAT_T score = 0;
    initialize_parameters();    

    // set verbosity
    if(CARP_FATAL <= verbosity && verbosity <= CARP_MAX){
      set_verbosity_level(verbosity);
    }
    else{
      carp(CARP_FATAL, "verbosity level must be between 0-100");
    }

    // set verbosity
    //set_verbosity_level(verbosity);
    
    if( charge < 1 || charge > 6){
      carp(CARP_FATAL, "peptide charge must be between 1 and 6.");
    }

    // check peptide sequence
    if(!valid_peptide_sequence(sequenceA)) {
      carp(CARP_ERROR, "%s not a valid sequence", sequenceA);
    }

    if(!valid_peptide_sequence(sequenceB) && strcmp(sequenceB, "NULL") != 0) {
      carp(CARP_ERROR, "%s not a valid sequence", sequenceB);
    }
    
    // parameters are now confirmed, can't be changed
    //parameters_confirmed();
    
    // create new ion series
    //LinkedIonSeries ion_series = LinkedIonSeries(positions, charge); 
    LinkedPeptide::linker_mass = atof(linker_mass_string);
    // a single peptide linked to itself
    if (strcmp(sequenceB, "NULL") == 0) {
      cout << "B is null" << endl; 
      sequenceB = NULL;
    }
    LinkedIonSeries ion_series;
    LinkedPeptide lp = LinkedPeptide(sequenceA, sequenceB, positionA, positionB, charge);
    //cout << lp << endl;
    ion_series.add_linked_ions(lp);
    carp(CARP_DEBUG, "number of ions: %d\n", ion_series.size());
    ion_series.print(); 

   // read ms2 file
    collection = new_spectrum_collection(ms2_file);
    spectrum = allocate_spectrum();
    //cout << "lp " << lp << endl; 
    // search for spectrum with correct scan number
    if(!get_spectrum_collection_spectrum(collection, scan_num, spectrum)){
      carp(CARP_ERROR, "failed to find spectrum with  scan_num: %d", scan_num);
      free_spectrum_collection(collection);
      free_spectrum(spectrum);
      exit(1);
    }
    Scorer xhhc_scorer = Scorer();
    if (print_spectrums) xhhc_scorer.set_print(true);
    else xhhc_scorer.set_print(false);
    score = xhhc_scorer.score_spectrum_vs_series(spectrum, ion_series);
    cout << "Xcorr score is: " << score << endl;  
   // free heap
   free_spectrum_collection(collection);
   free_spectrum(spectrum);
   free_parameters();
 }
 else{
   char* usage = parse_arguments_get_usage("xhhc-score-peptide-spectrum");
   result = parse_arguments_get_error(&error_message);
   fprintf(stderr, "Error in command line. Error # %d\n", result);
   fprintf(stderr, "%s\n", error_message);
   fprintf(stderr, "%s", usage);
   free(usage);
 }
 exit(0);
}

//creates three files for spectacle.pl: spectrums.out, theoretical.out, observed.out
