#include "crux_ion_series.h"

#define MAX_IONS 10000

void set_ion_series_is_predicted(ION_SERIES_T* ion_series, BOOLEAN_T is_predicted) {
  int offset = sizeof(char*) + //peptide
    sizeof(MODIFIED_AA_T*) + //modified_aa_seq
    sizeof(FLOAT_T) + //peptide_mass
    sizeof(int) + //charge
    sizeof(ION_CONSTRAINT_T*) + //constraint
    sizeof(ION_T*) * MAX_IONS + //ions
    sizeof(int); //num_ions

  char* ptr = (char*)ion_series + offset;

  *((BOOLEAN_T*)ptr) = is_predicted; 

}

BOOLEAN_T get_ion_series_is_predicted(ION_SERIES_T* ion_series) {
  int offset = sizeof(char*) + //peptide
    sizeof(MODIFIED_AA_T*) + //modified_aa_seq
    sizeof(FLOAT_T) + //peptide_mass
    sizeof(int) + //charge
    sizeof(ION_CONSTRAINT_T*) + //constraint
    sizeof(ION_T*) * MAX_IONS + //ions
    sizeof(int); //num_ions

  char* ptr = (char*)ion_series + offset;
  
  return *((BOOLEAN_T*)ptr);
}

MODIFIED_AA_T* get_ion_series_modified_aa_seq(ION_SERIES_T* ion_series) {
  int offset = sizeof(char*);
  char* ptr = (char*)ion_series + offset;

  return (MODIFIED_AA_T*)ptr;

}

BOOLEAN_T get_ion_constraint_use_neutral_losses(ION_CONSTRAINT_T* ion_constraint) {
  char* ptr = (char*)ion_constraint;

  return *((BOOLEAN_T*)ptr);
}


