#include "MPSM_MatchCollection.h"


MPSM_MatchCollection::MPSM_MatchCollection() {

}
  
BOOLEAN_T MPSM_MatchCollection::addMatch(MPSM_Match& match) {

//Add match if it is unique, i.e. we haven't seen this match yet.

  matches_.push_back(match);
  return TRUE;

}

MPSM_Match& MPSM_MatchCollection::getMatch(int idx) {
  return matches_[idx];
}
  
MPSM_Match& MPSM_MatchCollection::operator [] (int idx) {
  return getMatch(idx);
}
  



int MPSM_MatchCollection::numMatches() {
  return matches_.size();
}

void MPSM_MatchCollection::sortByScore(SCORER_TYPE_T match_mode) {
  //sort by score.
  for (int i=0;i<numMatches();i++) {
    matches_[i].getScore(match_mode);
  }
}
