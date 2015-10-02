/**
 * \file XLinkablePeptideIteratorTopN.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Iterator for xlinkable peptides
 *****************************************************************************/

#include "XLinkablePeptideIteratorTopN.h"
#include "XLinkablePeptideIterator.h"
#include "XLink.h"
#include "XLinkPeptide.h"
#include "XLinkScorer.h"
#include "XLinkDatabase.h"
#include "util/GlobalParams.h"
#include <iostream>


using namespace std;

/**
 * constructor that sets up the iterator
 */
XLinkablePeptideIteratorTopN::XLinkablePeptideIteratorTopN(
  Crux::Spectrum* spectrum, ///<spectrum
  FLOAT_T precursor_mass, ///< Mass of precursor
  FLOAT_T min_mass, ///< min mass of candidates
  FLOAT_T max_mass, ///< max mass of candidates
  int precursor_charge, ///< Charge of precursor  
  bool is_decoy ///< generate decoy candidates
  ) {

  carp(CARP_DEBUG, "XLinkablePeptideIteratorTopN: start()");

  scored_xlp_.clear();

  XLinkScorer scorer(spectrum, precursor_charge);
  top_n_ = GlobalParams::getXLinkTopN();
  //carp(CARP_DEBUG, "top_in:%i", top_n_);
  //carp(CARP_DEBUG, "precursor:%g", precursor_mass); 
  //carp(CARP_INFO, "min:%g", min_mass);
  //carp(CARP_INFO, "max:%g", max_mass);
  
  //find the begin and end iterators for the given mass range.
  vector<XLinkablePeptide>::iterator biter = XLinkDatabase::getXLinkableFlattenBegin(is_decoy,min_mass);
  vector<XLinkablePeptide>::iterator eiter = XLinkDatabase::getXLinkableFlattenEnd(is_decoy, max_mass);

  //carp(CARP_INFO, "biter:%d", biter-XLinkDatabase::getXLinkableFlattenBegin());
  //carp(CARP_INFO, "eiter:%d", eiter-XLinkDatabase::getXLinkableFlattenBegin());
  //carp(CARP_INFO, "eiter:%f", eiter->getMass(MONO));

  //check the range to see if we need to score.
  int range = eiter - biter;
  //carp(CARP_INFO, "range:%d", range);

  if (range > top_n_) {

    while(biter != eiter) {
      XLinkablePeptide& pep1 = *biter;
      FLOAT_T delta_mass = precursor_mass - pep1.getMass(MONO) - XLinkPeptide::getLinkerMass();
      FLOAT_T xcorr = scorer.scoreXLinkablePeptide(pep1, 0, delta_mass);
      pep1.setXCorr(0, xcorr);
      scored_xlp_.push_back(pep1);
      biter++;
    }
    carp(CARP_DEBUG, "number of xlinkable peptides scored:%d", scored_xlp_.size());
    sort(scored_xlp_.begin(), scored_xlp_.end(), compareXLinkableXCorr);
  } else {
    //carp(CARP_INFO, "No scoring needed!");
    scored_xlp_.insert(scored_xlp_.begin(), biter, eiter);
  }
  current_count_ = 0;
  queueNextPeptide();
  
  

}

/**
 * Destructor
 */
XLinkablePeptideIteratorTopN::~XLinkablePeptideIteratorTopN() {
}

/**
 * queues the next linkable peptide
 */
void XLinkablePeptideIteratorTopN::queueNextPeptide() {
  //carp(CARP_DEBUG, "queueNextPeptide()::start");

  if (current_count_ < scored_xlp_.size() && current_count_ < top_n_) {
    //current_ = scored_xlp_[current_count_];
    current_count_++;
    has_next_ = true;
  } else {
    has_next_ = false;
  }
  //carp(CARP_DEBUG, "queueNextPeptide:%i", has_next_?1:0);
}

/**
 *\returns whether there is another linkable peptide
 */
bool XLinkablePeptideIteratorTopN::hasNext() {
  //carp(CARP_DEBUG, "XLinkablePeptideIteratorTopN::hasNext(): %i", has_next_?1:0);
  return has_next_;
}

/**
 *\returns the next peptide
 */
XLinkablePeptide& XLinkablePeptideIteratorTopN::next() {
  //carp(CARP_INFO, "XLinkablePeptideIteratorTopN::next()");
  if (!has_next_) {
    carp(CARP_FATAL, "next called on empty iterator!");
  }

  XLinkablePeptide& ans = scored_xlp_[current_count_-1];
  //carp(CARP_INFO, "next peptide:%s %g", ans.getSequence(), ans.getXCorr());
  queueNextPeptide();
  //carp(CARP_INFO, "XLinkablePeptideIteratorTopN: returning reference");
  return ans;
}

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
