
//#include "hhc.h"
#include "xhhc_ion_series.h"
//#include "xhhc.h"

#define NUM_ARGUMENTS 7
#define NUM_OPTIONS 0
using namespace std;

int main(int argc, char** argv) {

  char* peptideA = NULL;
  char* peptideB = NULL;
  int posA = 0;
  int posB = 0;
  int linker_mass = 0;
  int charge = 1; 

 parse_arguments_set_req(
	"peptide A",
 	"first peptide sequence",
	(void *) &peptideA,
	STRING_ARG);

 parse_arguments_set_req(
	"peptide B",
 	"second peptide sequence",
	(void *) &peptideB,
	STRING_ARG);

 parse_arguments_set_req(
	"position A",
 	"zero based index of linker on peptide A [0, length(A)-1]",
	(void *) &posA,
	INT_ARG);

 parse_arguments_set_req(
	"position B",
 	"zero based index of linker on peptide B [0, length(B)-1]",
	(void *) &posB,
	INT_ARG);

 parse_arguments_set_req(
	"charge",
 	"linked peptide charge",
	(void *) &charge,
	INT_ARG);

 parse_arguments_set_req(
	"linker-mass",
 	"mass of the link between A and B",
	(void *) &linker_mass,
	INT_ARG);

  set_verbosity_level(CARP_INFO);
  if (!parse_arguments(argc, argv, 0)) {
   char* error_message;
   char* usage = parse_arguments_get_usage("xhhc-predict-peptide-ions");
   int result = parse_arguments_get_error(&error_message);
   fprintf(stderr, "Error in command line. Error # %d\n", result);
   fprintf(stderr, "%s\n", error_message);
   fprintf(stderr, "%s", usage);
   free(usage);
   exit(1);
 }

  LinkedPeptide::linker_mass = linker_mass; 
  LinkedPeptide linked_peptide;

  if (string(peptideB) == "NULL")
    linked_peptide = LinkedPeptide(peptideA, NULL, posA, posB, charge);
  else  
    linked_peptide = LinkedPeptide( peptideA, peptideB, posA, posB, charge);  

  cout << "precursor: " << linked_peptide << endl;
  LinkedIonSeries ion_series;
  ion_series.add_linked_ions(linked_peptide);
  ion_series.print(); 
  return 0;
}
