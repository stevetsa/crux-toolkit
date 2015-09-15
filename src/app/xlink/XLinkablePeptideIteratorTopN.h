/**
 * \file XLinkablePeptideIterator.h
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Iterator for xlinkable peptides
 *****************************************************************************/

#ifndef XLINKABLEPEPTIDEITERATORTOPN_H
#define XLINKABLEPEPTIDEITERATORTOPN_H

#include "XLinkablePeptideIterator.h"
#include "XLinkablePeptide.h"
#include "XLinkBondMap.h"

#include <vector>

class XLinkablePeptideIteratorTopN: public XLinkablePeptideIterator {

 protected:


  std::vector<XLinkablePeptide> scored_xlp_; ///< sorted by highest XCorr score.
  int current_idx_; 
  int current_rank_;
  FLOAT_T current_xcorr_;
  int top_n_; ///<set by kojak-top-n
  bool has_next_; ///< is there a next candidate
  bool is_decoy_; ///< are we getting decoys

  /**
   * queues the next linkable peptide
   */
  void queueNextPeptide(); 

 public:

  /**
   * constructor that sets up the iterator
   */
  XLinkablePeptideIteratorTopN(
    Crux::Spectrum* spectrum,
    FLOAT_T precursor_mass,
    FLOAT_T min_mass, ///< min mass of candidates
    FLOAT_T max_mass, ///< max mass of candidates
    int precursor_charge, ///< Charge of precursor
    Database* database, ///<peptide index
    PEPTIDE_MOD_T** peptide_mods, ///< current peptide mod
    int num_peptide_mods, 
    bool is_decoy, ///< generate decoy candidates
    XLinkBondMap& bondmap ///< map of valid links
    );

  /**
   * Destructor
   */
  virtual ~XLinkablePeptideIteratorTopN();

  /**
   * \returns whether there is another linkable peptide
   */
  bool hasNext();

  /**
   *\returns the next linkable peptide
   */
  XLinkablePeptide next();

};


#endif

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
