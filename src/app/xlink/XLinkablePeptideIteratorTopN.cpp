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
      onelink.getXCorr();
      scored_xlp_.push_back(onelink);
    }
  }
  carp(CARP_DEBUG, "number of xlinkable peptides scored:%d", scored_xlp_.size());

  sort(scored_xlp_.begin(), scored_xlp_.end(), compareXLinkableXCorr);
  
  for (size_t idx=0;idx < scored_xlp_.size() ;idx++) {
    //carp(CARP_DEBUG, "%i:%s xcorr:%g", idx, scored_xlp_[idx].getSequence(), scored_xlp_[idx].getXCorr());
  }
  current_idx_ = -1;
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
  if (current_idx_ == -1) {
    carp(CARP_DEBUG, "first call");
    if (scored_xlp_.empty()) {
      carp(CARP_DEBUG, "empty scored peptides");
      has_next_ = false;
    } else {
      carp(CARP_DEBUG, "Starting at first peptide");
      current_ = scored_xlp_.at(0);
      current_xcorr_ = current_.getXCorr();
      current_rank_ = 1;
      has_next_ = true;
      current_idx_ = 0; 
    }
  } else if (current_idx_ < (scored_xlp_.size()-1)) {
    size_t next_idx = current_idx_ + 1;
    size_t next_rank;
    FLOAT_T next_xcorr = scored_xlp_.at(next_idx).getXCorr();
    if (next_xcorr == current_xcorr_) {
      //rank stays the same.
      next_rank = current_rank_;
    } else {
      //rank is the idx.
      next_rank = next_idx+1;
    }
    if (next_rank <= top_n_) {
      //Success!
      has_next_ = true;
      current_idx_ = next_idx;
      current_xcorr_ = next_xcorr;
      current_rank_ = next_rank;
      current_ = scored_xlp_.at(current_idx_);
    } else {
      has_next_ = false;
    }
  } else {
    has_next_ = false;
  }
  carp(CARP_DEBUG, "queueNextPeptide:%i", has_next_?1:0);
}

/**
 *\returns whether there is another linkable peptide
 */
bool XLinkablePeptideIteratorTopN::hasNext() {
  carp(CARP_DEBUG, "XLinkablePeptideIteratorTopN::hasNext(): %i", has_next_?1:0);
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
  //carp(CARP_DEBUG, "next peptide:%s %g", ans.getSequence(), ans.getXCorr());
  queueNextPeptide();
  return ans;
}

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
