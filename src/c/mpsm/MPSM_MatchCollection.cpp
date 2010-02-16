#include "MPSM_MatchCollection.h"



MPSM_MatchCollection::MPSM_MatchCollection() {

}

MPSM_MatchCollection::MPSM_MatchCollection(MATCH_COLLECTION_T* spsm_matches) {
  MATCH_ITERATOR_T *match_iter =
    new_match_iterator(spsm_matches, XCORR, FALSE);

  while(match_iterator_has_next(match_iter)) {
    MATCH_T* match = match_iterator_next(match_iter);
    MPSM_Match mpsm_match(match);
    addMatch(mpsm_match);
  }

  free_match_iterator(match_iter);
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
