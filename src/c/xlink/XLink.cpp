#include "XLink.h"


#include <sstream>

using namespace std;

namespace XLink {

set<PEPTIDE_T*> allocated_peptides_;

void addAllocatedPeptide(PEPTIDE_T* peptide) {

  allocated_peptides_.insert(peptide);
}

void deleteAllocatedPeptides() {
  for (set<PEPTIDE_T*>::iterator iter =
    allocated_peptides_.begin();
    iter != allocated_peptides_.end();
    ++iter) {
  
  free_peptide(*iter);

  }
  allocated_peptides_.clear();
}


void get_protein_ids_locations(PEPTIDE_T *peptide, 
  set<string>& protein_ids_locations) {

  PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator = 
    new_peptide_src_iterator(peptide);

  std::ostringstream protein_field_stream;

  if (peptide_src_iterator_has_next(peptide_src_iterator)) {
    while(peptide_src_iterator_has_next(peptide_src_iterator)){
      PEPTIDE_SRC_T* peptide_src = peptide_src_iterator_next(peptide_src_iterator);
      Protein* protein = get_peptide_src_parent_protein(peptide_src);
      char* protein_id = protein->getId();
      int peptide_loc = get_peptide_src_start_idx(peptide_src);
      std::ostringstream protein_loc_stream;
      protein_loc_stream << protein_id << "(" << peptide_loc << ")";
      free(protein_id);
      protein_ids_locations.insert(protein_loc_stream.str());
    }
  }
  free(peptide_src_iterator);
}

string get_protein_ids_locations(PEPTIDE_T* peptide) {
  set<string> protein_ids_locations;

  get_protein_ids_locations(peptide, protein_ids_locations);


  set<string>::iterator result_iter = protein_ids_locations.begin();

  string protein_field_string = *result_iter;

  while(++result_iter != protein_ids_locations.end()) {
    protein_field_string += "," + *result_iter;
  }

  return protein_field_string;

}





};
