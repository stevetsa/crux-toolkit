/**
 * \file XLinkablePeptideIterator.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Iterator for xlinkable peptides
 *****************************************************************************/

#include "XLinkablePeptideIterator.h"
#include "XLink.h"
#include "LinearPeptide.h"
#include "util/GlobalParams.h"
#include <iostream>

using namespace std;

std::vector<XLinkablePeptide> XLinkablePeptideIterator::linkable_peptides_vec_;
std::vector<XLinkablePeptide> XLinkablePeptideIterator::decoy_linkable_peptides_vec_;


void XLinkablePeptideIterator::generateAllLinkablePeptides(
  vector<XLinkablePeptide>& xlp,
  Database* database,
  PEPTIDE_MOD_T** peptide_mods,
  int num_peptide_mods, 
  bool decoy) {

  xlp.empty();

  int max_mod_xlink = GlobalParams::getMaxXLinkMods();

  for (int mod_idx=0;mod_idx < num_peptide_mods; mod_idx++) {
    PEPTIDE_MOD_T* peptide_mod = peptide_mods[mod_idx];
    ModifiedPeptidesIterator* peptide_iterator =
      new ModifiedPeptidesIterator(
	  			 GlobalParams::getMinMass(),
		  		 GlobalParams::getMaxMass(),
			  	 peptide_mod,
				 false,
				 database,
				 1);

    while (peptide_iterator->hasNext()) {
      Crux::Peptide* peptide = peptide_iterator->next();
      if (peptide->countModifiedAAs() <= max_mod_xlink) {
        XLinkablePeptide::findLinkSites(peptide, bondmap_, link_sites_);
        if (link_sites_.size() > 0) {
          has_next_ = true;
          XLinkablePeptide current = XLinkablePeptide(peptide, link_sites_);
          current.setDecoy(decoy);
          xlp.push_back(current);  
	}       
      }
    }
    delete peptide_iterator;
  }

  sort(xlp.begin(), xlp.end(), compareXLinkablePeptideMass);
}

/**
 * constructor that sets up the iterator
 */
XLinkablePeptideIterator::XLinkablePeptideIterator(
    double min_mass, ///< min mass of candidates
    double max_mass, ///< max mass of candidates
    Database* database, ///< protein database
    PEPTIDE_MOD_T** peptide_mods, ///<current peptide mod
    int num_peptide_mods,
    bool is_decoy, ///< generate decoy candidates
    XLinkBondMap& bondmap ///< map of valid links
    ) {

  is_decoy_ = is_decoy;

  bondmap_ = bondmap;
  min_mass_ = min_mass;
  max_mass_ = max_mass;
  if (is_decoy_) {
    if (decoy_linkable_peptides_vec_.empty()) {
      generateAllLinkablePeptides(decoy_linkable_peptides_vec_, database, 
        peptide_mods, num_peptide_mods, is_decoy_);
    }
    iter_ = lower_bound(decoy_linkable_peptides_vec_.begin(), decoy_linkable_peptides_vec_.end(), 
      min_mass_, compareXLinkablePeptideMassToFLOAT);
    has_next_ = iter_ != decoy_linkable_peptides_vec_.end() &&
      iter_->getMass() <= max_mass_;
  } else {
    if (linkable_peptides_vec_.empty()) {
      generateAllLinkablePeptides(linkable_peptides_vec_, database, peptide_mods, num_peptide_mods, is_decoy_);
    }
    iter_ = lower_bound(linkable_peptides_vec_.begin(), linkable_peptides_vec_.end(), 
      min_mass_, compareXLinkablePeptideMassToFLOAT);
    has_next_ = iter_ != linkable_peptides_vec_.end() && iter_->getMass() <= max_mass_;
  }
}

/**
 * Destructor
 */
XLinkablePeptideIterator::~XLinkablePeptideIterator() {
  //delete peptide_iterator_;
}

/**
 * queues the next linkable peptide
 */
void XLinkablePeptideIterator::queueNextPeptide() {
  
  iter_++;
  if (is_decoy_) {
    has_next_ = iter_ != decoy_linkable_peptides_vec_.end() && iter_->getMass() <= max_mass_;
  } else {
    has_next_ = iter_ != linkable_peptides_vec_.end() && iter_->getMass() <= max_mass_;
  }
}

/**
 *\returns whether there is another linkable peptide
 */
bool XLinkablePeptideIterator::hasNext() {

  return has_next_;
}

/**
 *\returns the next peptide
 */
XLinkablePeptide XLinkablePeptideIterator::next() {

  if (!has_next_) {
    carp(CARP_WARNING, "next called on empty iterator!");
  }
  
  XLinkablePeptide ans = *iter_;
  queueNextPeptide();
  return ans;
}

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
