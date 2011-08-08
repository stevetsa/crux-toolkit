#ifndef LINEARPEPTIDE_H_
#define LINEARPEPTIDE_H_

#include "objects.h"
#include "utils.h"

#include "XLinkMatch.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <vector>

class LinearPeptide : public XLinkMatch {
 protected:

  PEPTIDE_T* peptide_;
  char* sequence_;
  BOOLEAN_T is_decoy_;
 public:
  
  LinearPeptide();
  LinearPeptide(char* sequence_);
  LinearPeptide(PEPTIDE_T* peptide);

  virtual ~LinearPeptide();

  static void addCandidates(
    FLOAT_T min_mass,
    FLOAT_T max_mass,
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
