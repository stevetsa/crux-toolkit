#include "XLink.h"

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

};
