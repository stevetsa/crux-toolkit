#ifndef XLINKPEPTIDE_H_
#define XLINKPEPTIDE_H_

#include "objects.h"
#include "utils.h"

#include "XLinkMatch.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <set>
#include <vector>

class XLinkPeptide : public XLinkMatch {
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

  bool isInter();
  bool isIntra();


  static void setLinkerMass(FLOAT_T linker_mass);
  static FLOAT_T getLinkerMass();
  static void addCandidates(
    FLOAT_T min_mass, 
    FLOAT_T max_mass,
    XLinkBondMap& bondmap, 
    Index* index, 
    Database* database,
    PEPTIDE_MOD_T** peptide_mods,
    int num_peptide_mods,
    XLinkMatchCollection& candidates);

  static void addLinkablePeptides(
    double min_mass, double max_mass,
    Index* index, Database* database,
    PEPTIDE_MOD_T* peptide_mod, BOOLEAN_T is_decoy, 
    XLinkBondMap& bondmap, 
    std::vector<XLinkablePeptide>& linkable_peptides);

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
