#ifndef XLINKMATCHCOLLECTION_H_
#define XLINKMATCHCOLLECTION_H_

/* Crux Includes */
#include "objects.h"
#include "MatchCollection.h"
#include "Index.h"
#include "Database.h"
#include "modifications.h"
#include "SpectrumZState.h"

/* XLink Includes */
#include "XLink.h"
#include "XLinkMatch.h"

class XLinkMatchCollection : public MatchCollection {
 protected:

  bool include_linear_peptides;
  bool include_self_loops;

  SpectrumZState zstate_;
  int scan_;
  FLOAT_T precursor_mz_;
  Spectrum* spectrum_;

  void addCandidates(
    FLOAT_T min_mass,
    FLOAT_T max_mass,
    XLinkBondMap& bondmap,
    Index* index,
    Database* database,
    PEPTIDE_MOD_T** peptide_mods,
    int num_peptide_mods);


 public:
  XLinkMatchCollection();
  XLinkMatchCollection(XLinkMatchCollection& vector);

  XLinkMatchCollection(
    XLinkBondMap& bondmap,
    PEPTIDE_MOD_T** peptide_mods,
    int num_peptide_mods,
    Index* index,
    Database* database);


  XLinkMatchCollection(FLOAT_T precursor_mz,
                       SpectrumZState& zstate,
		       XLinkBondMap& bondmap,
		       Index* index,
		       Database* database,
		       PEPTIDE_MOD_T** peptide_mods,
		       int num_peptide_mods,
		       bool is_decoy=FALSE);

  virtual ~XLinkMatchCollection();
  
  void add(XLinkMatch* candidate);


  void shuffle();
  void shuffle(XLinkMatchCollection& decoy_vector);

  void scoreSpectrum(Spectrum* spectrum);
  void sortByXCorr();
  void sortBySP();
  void setRanks();
  void fitWeibull(FLOAT_T& shift, 
		  FLOAT_T& eta, 
		  FLOAT_T& beta, 
		  FLOAT_T& corr);

  void setScan(unsigned int scan);
  unsigned int getScan();
  int getCharge();
  
  FLOAT_T getPrecursorMZ();
  FLOAT_T getSpectrumNeutralMass();

  unsigned int size() {return 0;} //PLACEHOLDER


};

#endif
