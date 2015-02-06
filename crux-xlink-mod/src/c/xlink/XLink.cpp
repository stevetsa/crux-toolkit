#include "XLink.h"


#include <sstream>

using namespace std;

namespace XLink {

set<Peptide*> allocated_peptides_;

void addAllocatedPeptide(Peptide* peptide) {

  allocated_peptides_.insert(peptide);
}

void deleteAllocatedPeptides() {
  for (set<Peptide*>::iterator iter =
    allocated_peptides_.begin();
    iter != allocated_peptides_.end();
    ++iter) {
  
  delete *iter;

  }
  allocated_peptides_.clear();
}


void get_protein_ids_locations(Peptide *peptide, 
  set<string>& protein_ids_locations) {

  std::ostringstream protein_field_stream;

  for (PeptideSrcIterator peptide_src_iterator =
    peptide->getPeptideSrcBegin();
    peptide_src_iterator != peptide->getPeptideSrcEnd();
    ++peptide_src_iterator) {
    
    PeptideSrc* peptide_src = *peptide_src_iterator;
    Protein* protein = peptide_src->getParentProtein();
    char* protein_id = protein->getId();
    int peptide_loc = peptide_src->getStartIdx();
    std::ostringstream protein_loc_stream;
    protein_loc_stream << protein_id << "(" << peptide_loc << ")";
    free(protein_id);
    protein_ids_locations.insert(protein_loc_stream.str());
  }

}

string get_protein_ids_locations(Peptide* peptide) {
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
