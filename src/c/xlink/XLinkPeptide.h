#ifndef XLINKPEPTIDE_H_
#define XLINKPEPTIDE_H_

#include "objects.h"
#include "utils.h"

#include "MatchCandidate.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <set>
#include <vector>

class XLinkPeptide : public MatchCandidate {
 protected:
  static FLOAT_T linker_mass_;
  static std::set<PEPTIDE_T*> allocated_peptides_;
  std::vector<XLinkablePeptide> linked_peptides_;
  std::vector<int> link_pos_idx_;
  
  bool mass_calculated_[NUMBER_MASS_TYPES];
  FLOAT_T mass_[NUMBER_MASS_TYPES];
  
  BOOLEAN_T is_decoy_;
  
  int getLinkPos(int peptide_idx);

 public:
  
  XLinkPeptide();
  XLinkPeptide(XLinkablePeptide& peptideA, 
	       XLinkablePeptide& peptideB, 
	       int posA, int posB);
  XLinkPeptide(char* peptideA,
	       char* peptideB,
	       int posA, int posB);

  virtual ~XLinkPeptide();

  static void setLinkerMass(FLOAT_T linker_mass);
  static FLOAT_T getLinkerMass();
  static void addCandidates(FLOAT_T precursor_mz, int charge, 
			    XLinkBondMap& bondmap, 
			    INDEX_T* index, DATABASE_T* database,
			    PEPTIDE_MOD_T** peptide_modes,
			    int num_peptide_mods,
			    MatchCandidateVector& candidates,
			    BOOLEAN_T use_decoy_window=FALSE);

  static void addLinkablePeptides(
    double min_mass, double max_mass,
    INDEX_T* index, DATABASE_T* database,
    PEPTIDE_MOD_T* peptide_mod, BOOLEAN_T is_decoy, 
    XLinkBondMap& bondmap, 
    std::vector<XLinkablePeptide>& linkable_peptides);

  virtual MATCHCANDIDATE_TYPE_T getCandidateType();
  virtual std::string getSequenceString();
  virtual FLOAT_T calcMass(MASS_TYPE_T mass_type);

  virtual MatchCandidate* shuffle();

  virtual void predictIons(ION_SERIES_T* ion_series, int charge);
  std::string getIonSequence(ION_T* ion);
  virtual PEPTIDE_T* getPeptide(int peptide_idx);
};



#endif
