#ifndef SELFLOOPPEPTIDE_H_
#define SELFLOOPPEPTIDE_H_

#include "objects.h"
#include "utils.h"

#include "MatchCandidate.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <vector>

class SelfLoopPeptide : public MatchCandidate {
 protected:

  XLinkablePeptide linked_peptide_;
  std::vector<int> link_pos_idx_;
  
  BOOLEAN_T is_decoy_;
 public:
  
  int getLinkPos(int link_idx);

  SelfLoopPeptide();
  SelfLoopPeptide(XLinkablePeptide& peptide,
		  int posA, int posB);

  SelfLoopPeptide(
    char* peptide,
    int posA,
    int posB
  );

  virtual ~SelfLoopPeptide();

  static void addCandidates(FLOAT_T precursor_mz, int charge, 
			    XLinkBondMap& bondmap, 
			    INDEX_T* index, DATABASE_T* database,
			    PEPTIDE_MOD_T** peptide_mods,
			    int num_peptide_mods,
			    MatchCandidateVector& candidates,
			    BOOLEAN_T use_decoy_window=FALSE);

  virtual MATCHCANDIDATE_TYPE_T getCandidateType();
  virtual std::string getSequenceString();
  virtual FLOAT_T calcMass(MASS_TYPE_T mass_type);
  virtual MatchCandidate* shuffle();
  virtual void predictIons(ION_SERIES_T* ion_series, int charge);
  std::string getIonSequence(ION_T* ion);
  virtual PEPTIDE_T* getPeptide(int peptide_idx);
};



#endif
