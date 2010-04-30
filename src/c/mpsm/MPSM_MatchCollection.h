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
  void calcZScores();

  friend std::ostream& operator<<(std::ostream& os, MPSM_MatchCollection& collection_obj);


};

#endif
