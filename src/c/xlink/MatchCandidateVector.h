#ifndef MATCHCANDIDATEVECTOR_H_
#define MATCHCANDIDATEVECTOR_H_

#include "XLink.h"
#include "MatchCandidate.h"
#include "SpectrumZState.h"

#include "Index.h"
#include "Database.h"
#include "modifications.h"

#include <vector>

class MatchCandidateVector : public std::vector<MatchCandidate*> {
 protected:
  BOOLEAN_T include_linear_peptides;
  BOOLEAN_T include_self_loops;

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
  MatchCandidateVector();
  MatchCandidateVector(MatchCandidateVector& vector);

  MatchCandidateVector(
    XLinkBondMap& bondmap,
    PEPTIDE_MOD_T** peptide_mods,
    int num_peptide_mods,
    Index* index,
    Database* database);


  MatchCandidateVector(FLOAT_T precursor_mz,
                       SpectrumZState& zstate,
		       XLinkBondMap& bondmap,
		       Index* index,
		       Database* database,
		       PEPTIDE_MOD_T** peptide_mods,
		       int num_peptide_mods,
		       BOOLEAN_T is_decoy=FALSE);

  virtual ~MatchCandidateVector();
  
  void add(MatchCandidate* candidate);


  void shuffle();
  void shuffle(MatchCandidateVector& decoy_vector);

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




};

#endif
