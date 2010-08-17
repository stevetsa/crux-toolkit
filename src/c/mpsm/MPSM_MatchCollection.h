#ifndef MPSM_MATCHCOLLECTION_H
#define MPSM_MATCHCOLLECTION_H

#include <set>


#include "MPSM_Match.h"

#include "match.h"
#include "match_collection.h"



class MPSM_MatchCollection {

protected:
  MATCH_COLLECTION_T* spsm_matches_;

  std::vector<MPSM_Match> matches_;

  double zscore_mean_;
  double zscore_std_;
  double xcorr_2_;
  bool sorted_;

public:

  MPSM_MatchCollection();
  MPSM_MatchCollection(MATCH_COLLECTION_T* spsm_matches);
  virtual ~MPSM_MatchCollection();
  
  void free();

  
  bool addMatch(MPSM_Match& match);

  MPSM_Match& getMatch(int idx);
  
  MPSM_Match& operator [] (int idx);

  int numMatches();
  
  void sortByScore(SCORER_TYPE_T match_mode);

  void calcDeltaCN();
  double calcDeltaCNMatch(double xcorr);
  void calcZParameters(double& mean, double& std);
  void setZParameters(double mean, double std);
  double calcZScore(double xcorr);

  friend std::ostream& operator<<(std::ostream& os, MPSM_MatchCollection& collection_obj);


};

#endif
