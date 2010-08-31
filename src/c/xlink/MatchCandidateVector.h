#ifndef MATCHCANDIDATEVECTOR_H_
#define MATCHCANDIDATEVECTOR_H_

#include "xlink.h"
#include "MatchCandidate.h"

#include "index.h"
#include "database.h"
#include "modifications.h"

#include <vector>

class MatchCandidateVector : public std::vector<MatchCandidate*> {
 protected:
  BOOLEAN_T include_linear_peptides;
  BOOLEAN_T include_self_loops;
  int charge_;
  int scan_;
  FLOAT_T precursor_mz_;

 public:
  MatchCandidateVector();
  MatchCandidateVector(MatchCandidateVector& vector);
  MatchCandidateVector(FLOAT_T precursor_mz,
		       int charge,
		       XLinkBondMap& bondmap,
		       INDEX_T* index,
		       DATABASE_T* database,
		       PEPTIDE_MOD_T** peptide_mods,
		       int num_peptide_mods,
		       BOOLEAN_T is_decoy=FALSE);

  virtual ~MatchCandidateVector();
  
  void add(MatchCandidate* candidate);


  void shuffle();
  void shuffle(MatchCandidateVector& decoy_vector);

  void scoreSpectrum(SPECTRUM_T* spectrum);
  void sortByXCorr();
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
