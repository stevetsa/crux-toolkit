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
  Database* database, ///< protein database
  PEPTIDE_MOD_T** peptide_mods, ///<current peptide mod
  int num_peptide_mods,
  bool is_decoy, ///< generate decoy candidates
  XLinkBondMap& bondmap ///< map of valid links
  ) {

  carp(CARP_DEBUG, "XLinkablePeptideIteratorTopN: start()");

  //scored_xlp_ = priority_queue<XLinkablePeptide,compareXLinkableXCorr>();

  XLinkScorer scorer(spectrum, precursor_charge);
  top_n_ = get_int_parameter("xlink-top-n");
  carp(CARP_DEBUG, "top_in:%i", top_n_);
  carp(CARP_DEBUG, "precursor:%g", precursor_mass); 
  carp(CARP_DEBUG, "min:%g", min_mass);
  carp(CARP_DEBUG, "max:%g", max_mass);

  vector<XLinkablePeptide>::iterator biter = XLinkDatabase::getXLinkableFlattenBegin(min_mass);
  vector<XLinkablePeptide>::iterator eiter = XLinkDatabase::getXLinkableFlattenEnd();

  while(biter != eiter && biter->getMass(MONO) <= max_mass) {
    XLinkablePeptide& pep1 = *biter;
    FLOAT_T delta_mass = precursor_mass - pep1.getMass(MONO) - XLinkPeptide::getLinkerMass();
    FLOAT_T xcorr = scorer.scoreXLinkablePeptide(pep1, 0, delta_mass);
    pep1.setXCorr(0, xcorr);
    scored_xlp_.push_back(pep1);
    biter++;
  }
  //carp(CARP_INFO, "number of xlinkable peptides scored:%d", scored_xlp_.size());
  if (scored_xlp_.size() > top_n_) {
    //dont sort if we are getting all of the candidates
    sort(scored_xlp_.begin(), scored_xlp_.end(), compareXLinkableXCorr);
  }
  //for (size_t idx=0;idx < scored_xlp_.size() ;idx++) {
    //carp(CARP_DEBUG, "%i:%s xcorr:%g", idx, scored_xlp_[idx].getSequence(), scored_xlp_[idx].getXCorr());
  //}
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
