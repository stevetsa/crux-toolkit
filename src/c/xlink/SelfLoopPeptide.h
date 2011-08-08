#ifndef SELFLOOPPEPTIDE_H_
#define SELFLOOPPEPTIDE_H_

#include "objects.h"
#include "utils.h"

#include "XLinkMatch.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <vector>

class SelfLoopPeptide : public XLinkMatch {
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

  static void addCandidates(
    FLOAT_T min_mass,
    FLOAT_T max_mass,
    XLinkBondMap& bondmap, 
    Index* index, 
    Database* database,
    PEPTIDE_MOD_T** peptide_mods,
    int num_peptide_mods,
    XLinkMatchCollection& candidates);

  virtual XLINKMATCH_TYPE_T getCandidateType();
  virtual std::string getSequenceString();
  virtual FLOAT_T calcMass(MASS_TYPE_T mass_type);
  virtual XLinkMatch* shuffle();
  virtual void predictIons(IonSeries* ion_series, int charge);
  std::string getIonSequence(Ion* ion);
  virtual PEPTIDE_T* getPeptide(int peptide_idx);

  virtual int getNumMissedCleavages();

  virtual bool isModified();

};



#endif
