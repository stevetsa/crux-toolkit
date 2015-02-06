#ifndef XLINKABLEPEPTIDEITERATOR_H
#define XLINKABLEPEPTIDEITERATOR_H

#include "ModifiedPeptidesIterator.h"
#include "XLinkablePeptide.h"
#include "XLinkBondMap.h"

class XLinkablePeptideIterator {

 protected:
  ModifiedPeptidesIterator* peptide_iterator_;
  XLinkablePeptide current_;
  XLinkBondMap bondmap_;
  bool has_next_;
  bool is_decoy_;
  std::vector<int> link_sites_;

  void queueNextPeptide();

 public:

  XLinkablePeptideIterator(
    double min_mass,
    double max_mass,
    Index* index,
    Database* database,
    PEPTIDE_MOD_T* peptide_mod,
    bool is_decoy,
    XLinkBondMap& bondmap);

  virtual ~XLinkablePeptideIterator();

  bool hasNext();
  XLinkablePeptide next();

};


#endif
