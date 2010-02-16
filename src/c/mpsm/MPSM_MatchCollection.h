#ifndef MPSM_MATCHCOLLECTION_H
#define MPSM_MATCHCOLLECTION_H

#include <vector>


#include "MPSM_Match.h"

#include "match.h"
#include "match_collection.h"



class MPSM_MatchCollection {

protected:

  std::vector<MPSM_Match> matches_;
public:

  MPSM_MatchCollection();
  MPSM_MatchCollection(MATCH_COLLECTION_T* spsm_matches);

  
  BOOLEAN_T addMatch(MPSM_Match& match);
  MPSM_Match& getMatch(int idx);
  
  MPSM_Match& operator [] (int idx);
  

  int numMatches();

  void sortByScore(SCORER_TYPE_T match_mode);

  


};

#endif
