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
  XLinkablePeptideIterator xlp_iterator(min_mass, max_mass, 
    database, peptide_mods, num_peptide_mods, is_decoy, bondmap);

  while(xlp_iterator.hasNext()) {
    XLinkablePeptide pep1 = xlp_iterator.next();
    FLOAT_T delta_mass = precursor_mass - pep1.getMass(MONO) - XLinkPeptide::getLinkerMass();
    for (unsigned int link1_idx=0;link1_idx < pep1.numLinkSites(); link1_idx++) {

      FLOAT_T xcorr = scorer.scoreXLinkablePeptide(pep1, link1_idx, delta_mass);
      XLinkablePeptide onelink(pep1);
      onelink.clearSites();
      onelink.addLinkSite(pep1.getLinkSite(link1_idx));
      onelink.setXCorr(0, xcorr);
      //carp(CARP_DEBUG, "%s %g", onelink.getSequence(), xcorr);
      //onelink.getXCorr();
      scored_xlp_.push(onelink);
      
    }
  }
  carp(CARP_DEBUG, "number of xlinkable peptides scored:%d", scored_xlp_.size());

  //sort(scored_xlp_.begin(), scored_xlp_.end(), compareXLinkableXCorr);
  
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
  carp(CARP_DEBUG, "queueNextPeptide()::start");

  if (!scored_xlp_.empty() && current_count_ <= top_n_) {
    current_count_++;
    current_ = scored_xlp_.top();
    scored_xlp_.pop();
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
XLinkablePeptide XLinkablePeptideIteratorTopN::next() {
  carp(CARP_DEBUG, "XLinkablePeptideIteratorTopN::next()");
  if (!has_next_) {
    carp(CARP_FATAL, "next called on empty iterator!");
  }

  XLinkablePeptide ans = current_;
  //carp(CARP_INFO, "next peptide:%s %g", ans.getSequence(), ans.getXCorr());
  queueNextPeptide();
  return ans;
}

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
