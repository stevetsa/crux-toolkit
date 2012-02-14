#ifndef MPSM_MATCHCOLLECTION_H
#define MPSM_MATCHCOLLECTION_H

#include <set>


#include "MPSM_Match.h"

#include "Match.h"
#include "MatchCollection.h"



class MPSM_MatchCollection {

protected:
  MatchCollection* spsm_matches_;

  std::vector<MPSM_Match> matches_;

  std::set<vector<Match*>,CompareMPSM_MatchVisited> visited_;

  double zscore_mean_;
  double zscore_std_;
  double xcorr_2_;
  bool sorted_;
  SCORER_TYPE_T sort_mode_;
  ZStateIndex zstate_index_;
  bool zstate_valid_;

public:

  MPSM_MatchCollection();
  MPSM_MatchCollection(MatchCollection* spsm_matches);
  virtual ~MPSM_MatchCollection();
  
  void free();
  void invalidate();
  ZStateIndex& getZStateIndex();  
  //ChargeIndex& getChargeIndex();
  
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

  double getSpectrumRTime();
  double getPredictedRTime(MPSM_Match& match);

  void calcRanks(SCORER_TYPE_T scorer_type);

  bool visited(MPSM_Match& match);


  friend std::ostream& operator<<(std::ostream& os, MPSM_MatchCollection& collection_obj);


};

#endif
