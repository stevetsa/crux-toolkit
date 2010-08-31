#include "XLink.h"

using namespace std;


void getLinkablePeptides(MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator, 
			 XLinkBondMap& bondmap, 
			 std::vector<XLinkablePeptide>& peptides) {

  //Iterator through all of the peptides, and collect the peptides
  //that have at least one link site available.
  while(modified_peptides_iterator_has_next(peptide_iterator)) {
    PEPTIDE_T* current_peptide = modified_peptides_iterator_next(peptide_iterator);
    //try creating a XLinkablePeptide.
    XLinkablePeptide xlink_peptide(current_peptide, bondmap);
    if (xlink_peptide.isLinkable()) {
      peptides.push_back(xlink_peptide);
    } else {
      free_peptide(current_peptide);
    }
  }
}

void getLinkablePeptides(FLOAT_T precursor_mz,
			 int charge,
			 PEPTIDE_MOD_T* peptide_mod,
			 BOOLEAN_T is_decoy,
			 INDEX_T* index,
			 DATABASE_T* database,
			 XLinkBondMap& bondmap,
			 std::vector<XLinkablePeptide>& peptides) {

  //create modified peptide iterator.
  MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator =
      new_modified_peptides_iterator_from_mz(precursor_mz,
                                             charge,
                                             peptide_mod, 
                                             is_decoy,
                                             index,
                                             database);
  
  getLinkablePeptides(peptide_iterator, bondmap, peptides);
  

  free_modified_peptides_iterator(peptide_iterator);
}
