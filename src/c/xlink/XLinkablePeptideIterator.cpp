#include "XLinkablePeptideIterator.h"
#include "XLink.h"
#include <iostream>

using namespace std;

XLinkablePeptideIterator::XLinkablePeptideIterator(
    double min_mass,
    double max_mass,
    Index* index,
    Database* database,
    PEPTIDE_MOD_T* peptide_mod,
    bool is_decoy,
    XLinkBondMap& bondmap) {

  is_decoy_ = is_decoy;

  bondmap_ = bondmap;

  //cerr << "min mass:"<<min_mass<<endl;
  //cerr << "max mass:"<<max_mass<<endl;

  peptide_iterator_ =     
    new ModifiedPeptidesIterator(
      min_mass, 
      max_mass,
      peptide_mod, 
      is_decoy,
      index, 
      database);
/*
  while (peptide_iterator_->hasNext()) {
    Peptide* peptide = peptide_iterator_->next();
    cerr << "peptide is:"<<peptide->getSequence()<<" "<<peptide->getPeptideMass()<<endl;
  }

  delete peptide_iterator_;
  peptide_iterator_ =     
    new ModifiedPeptidesIterator(
      min_mass, 
      max_mass,
      peptide_mod, 
      is_decoy,
      index, 
      database);

  */

  queueNextPeptide();

  

}

XLinkablePeptideIterator::~XLinkablePeptideIterator() {
  delete peptide_iterator_;

}


void XLinkablePeptideIterator::queueNextPeptide() {

  has_next_ = false;
  int max_mod_xlink = get_int_parameter("max-xlink-mods");
  while (peptide_iterator_->hasNext() && !has_next_) {

    Peptide* peptide = peptide_iterator_->next();
    //cerr << "peptide is:"<<peptide->getSequence()<<endl;  
    if (peptide->countModifiedAAs() <= max_mod_xlink) {
      XLinkablePeptide::findLinkSites(peptide, bondmap_, link_sites_);
      if (link_sites_.size() > 0) {
        has_next_ = true;
        current_ = XLinkablePeptide(peptide, link_sites_);
        current_.setDecoy(is_decoy_);
        XLink::addAllocatedPeptide(peptide);
      }
    } 
    if (!has_next_) {
      delete peptide;
    }
//    Peptide::free(peptide);
  }
}

bool XLinkablePeptideIterator::hasNext() {

  return has_next_;
}

XLinkablePeptide XLinkablePeptideIterator::next() {

  if (!has_next_) {
    carp(CARP_WARNING, "next called on empty iterator!");
  }

  XLinkablePeptide ans = current_;
  queueNextPeptide();
  return ans;
}
