#ifndef MPSM_MATCHCOLLECTION_H
#define MPSM_MATCHCOLLECTION_H

#include "MPSM_Match.h"



class MPSM_MatchCollection {

protected:

  std::vector<MPSM_Match> matches_;
public:

  MPSM_MatchCollection();
  
  BOOLEAN_T addMatch(MPSM_Match& match);
  MPSM_Match& getMatch(int idx);
  
  MPSM_Match& operator [] (int idx);
  

  int numMatches();

  void sortByScore(SCORER_TYPE_T match_mode);

  


};

#endif
